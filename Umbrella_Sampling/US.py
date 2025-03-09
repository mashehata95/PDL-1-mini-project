#!/usr/bin/env python

from openmm import *
from openmm.app import *
from openmm.unit import * 
from sys import stdout
from openmmplumed import PlumedForce
import argparse
import numpy as np
import time
import pytraj as pt

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--topology", help="topology file (prmtop)", required=True)
parser.add_argument("-c", "--coordinate", help="coordinate file (inpcrd)", required=True)
parser.add_argument("-s", "--structure", help="pdb file", required=True)
parser.add_argument("-nvt", "--nvt", help="Number of steps NVT equilibration", required=True)
parser.add_argument("-npt", "--npt", help="Number of steps NPT equilibration", required=True)
parser.add_argument("-l","--us" ,help="Number of steps for US production", required=True)
parser.add_argument("-t","--temperature" ,help="teperature in kelvin", required=True)
parser.add_argument("-gpuID", help="GPU ID", required=True)

args = parser.parse_args()

time_nvt = (float(args.nvt) * 1000) / 0.002
time_npt = (float(args.npt) * 1000) / 0.002
time_us = (float(args.us) * 1000) / 0.002

def EnergyMinimization():
    # Create the system
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds, removeCMMotion=True)
    
    # Define the integrator
    integrator = LangevinIntegrator(float(args.temperature)*kelvin, 1.0/picosecond, 0.002*picoseconds)

    # Set up the simulation platform
    platform = Platform.getPlatformByName('CPU')
    properties = {"Threads": "6"}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.getPositions())

    # Minimize the system's energy
    print("\t- Energy minimization started at {}".format(time.ctime()))
    simulation.minimizeEnergy()
    print("\t- Energy minimization finished at {}\n".format(time.ctime()))
    
    all_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    PDBFile.writeFile(prmtop.topology, all_positions, open("em.pdb", 'w'))
  
def Heating():
    pdb_em = PDBFile("em.pdb")
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds, removeCMMotion=True, rigidWater=False)  
    integrator = LangevinIntegrator(float(args.temperature)*kelvin, 1.0/picosecond, 0.002*picoseconds)
    
    # Applying position restraint on Ca atoms and ligand heavy atoms
    print("\t- Applying position restrain on protein Ca atoms and Glycans...")
    pt_system = pt.iterload("../protein.inpcrd", "../protein.prmtop")
    pt_topology = pt_system.top
    # Applying harmonic position restraint on Ca atoms and Carbohydrate molecules
    print("\t- Applying position restrain on protein Ca atoms and Carbohydrate molecules...")
    restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    restraint.addGlobalParameter('k', 50.0 * kilocalories_per_mole / nanometers**2)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    
    for atom in pdb.topology.atoms():
        if atom.name == "CA":
            restraint.addParticle(atom.index, pdb.positions[atom.index])
    
    system.addForce(restraint)
        
    # Equilibrate the system in the NVT ensemble
    platform = Platform.getPlatformByName('CUDA')
    properties = {"Precision": "mixed", 'DeviceIndex': str(args.gpuID)}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb_em.getPositions())
    simulation.context.setVelocitiesToTemperature(float(args.temperature)*kelvin)
    
    # Adding periodic condition
    if inpcrd.boxVectors is not None: 
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  
    # Running Heating
    simulation.reporters.append(StateDataReporter("heating.log", 500, step=True, kineticEnergy=True, potentialEnergy=True,
                                                  totalEnergy=True, temperature=True, volume=True))
  
    print("\t- Heating for {}ns started at {}".format(args.nvt, time.ctime()))
    simulation.step(time_nvt)
    print("\t- Heating finished at {}\n".format(time.ctime()))
    
    # Save rst7 and pdb file
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True, enforcePeriodicBox=True)
    with open("heating.rst7", 'w') as f:
        f.write(XmlSerializer.serialize(state))
    
    # Save PDB file
    all_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    PDBFile.writeFile(prmtop.topology, all_positions, open("heating.pdb", 'w'))

def Equilibration():
    # Creating system and defining integrator
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)
    system.addForce(MonteCarloBarostat(1.0 * bar, float(args.temperature)*kelvin, 25))
    integrator = LangevinIntegrator(float(args.temperature)*kelvin, 1.0/picosecond, 0.002*picoseconds)
    pdb = PDBFile("heating.pdb")
  
    # Applying position restraint on Ca atoms and ligand heavy atoms
    print("\t- Applying position restrain on protein Ca atoms and Glycans...")
    pt_system = pt.iterload("../protein.inpcrd", "../protein.prmtop")
    pt_topology = pt_system.top
    # Applying harmonic position restraint on Ca atoms and Carbohydrate molecules
    print("\t- Applying position restrain on protein Ca atoms and Carbohydrate molecules...")
    restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    restraint.addGlobalParameter('k', 50.0 * kilocalories_per_mole / nanometers**2)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    
    for atom in pdb.topology.atoms():
        if atom.name == "CA":
            restraint.addParticle(atom.index, pdb.positions[atom.index])
    
    system.addForce(restraint)

    # Equilibrate the system in the NPT ensemble
    platform = Platform.getPlatformByName('CUDA')
    properties = {"Precision": "mixed", 'DeviceIndex': str(args.gpuID)}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.getPositions())
    simulation.context.setVelocitiesToTemperature(float(args.temperature)*kelvin)
  
    # Adding periodic condition
    if inpcrd.boxVectors is not None: 
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
    # Reading system state and coordinates
    with open("heating.rst7", 'r') as f:
        state = XmlSerializer.deserialize(f.read())
    simulation.context.setState(state)
    simulation.context.setPositions(state.getPositions())
  
    # Running equilibration
    simulation.reporters.append(StateDataReporter("equilibration.log", 500, step=True, kineticEnergy=True, potentialEnergy=True,
                                                  totalEnergy=True, temperature=True, volume=True))
  
    print("\t- Equilibration for {}ns started at {}".format(args.npt, time.ctime()))
    simulation.step(time_npt)
    print("\t- Equilibration finished at {}\n".format(time.ctime()))
    
    # Save rst7 and pdb file
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True, enforcePeriodicBox=True)
    with open("equilibration.rst7", 'w') as f:
        f.write(XmlSerializer.serialize(state))
    
    # Save PDB file
    all_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    PDBFile.writeFile(prmtop.topology, all_positions, open("equilibration.pdb", 'w'))
  
def US():
    global time1
    pdb = PDBFile("equilibration.pdb")
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds, removeCMMotion=True)
    integrator = LangevinIntegrator(float(args.temperature)*kelvin, 1.0/picosecond, 0.002*picoseconds) # NVT ensemble
    system.addForce(MonteCarloBarostat(1.0 * bar, float(args.temperature)*kelvin, 25)) # NPT ensemble
    
    # Steered Molecular Dynamics (SMD) Setup
    plumed_script = open("plumed_file.dat", "r").read() # Reading the plumed input file
    force = PlumedForce(plumed_script)
    system.addForce(force)
    
    platform = Platform.getPlatformByName('CUDA')
    properties = {"Precision": "mixed", 'DeviceIndex': str(args.gpuID)}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    
    # Ensuring the existences of PBC 
    if inpcrd.boxVectors is not None:  # This will check if the system has periodic boundary condition
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  
    # Reading system state and coordinates
    with open("equilibration.rst7", 'r') as f:
        state = XmlSerializer.deserialize(f.read())
    simulation.context.setState(state)
    simulation.context.setPositions(state.getPositions())
    
    # Run US simulation
    print('\t- Running US...')
    simulation.reporters.append(DCDReporter('US.dcd', 500, enforcePeriodicBox=True))
    simulation.reporters.append(StateDataReporter("US.log", 500, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True))
    simulation.reporters.append(StateDataReporter(stdout, 500, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True))
    time1 = time.time()
    simulation.step(time_us)
    # Save each frame
    
    # Saving simulation state
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True, enforcePeriodicBox=True)
    with open("US.rst7", 'w') as f:
        f.write(XmlSerializer.serialize(state))

    # Saving final coordinates
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(), open("us_output.pdb", 'w'))

if __name__ == "__main__":
      
      print("- Parsing topology file\n")
      prmtop = AmberPrmtopFile(args.topology)
      
      print("- Parsing coordinate file\n")
      inpcrd = AmberInpcrdFile(args.inpcrd)

      print("- Parsing pdb file\n")
      pdb = PDBFile(args.structure)
      
      print(f"- Running energy minimization for {i}\n")
      EnergyMinimization()
      
      print(f"- Running heating for {i}\n")
      Heating()
      
      print(f"- Running equilibration for {i}\n")
      Equilibration()
      
      print(f"- Running US for {i}\n")
      US()
      
      print('US Simulation is Completed!')
      print("US {} took {} minutes".format(time_us, round((time.time() - time1) / 60, 2)))

