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
parser.add_argument("-p","--topology",help="Amber topology file",required=True)
parser.add_argument("-c","--coordinates",help="Amber coordinates file",required=True)
parser.add_argument("-s","--structure",help="PDB file",required=True)
parser.add_argument("-npt","--npt",help="Number of steps NPT equilibration",required=True)
parser.add_argument("-smd",help="Number of steps for SMD production",required=True)
parser.add_argument("-t","--temperature", help="temperature in kelvin",required=True)
parser.add_argument("-gpuID",help="GPU ID",required=True)

args = parser.parse_args()

time_npt = (float(args.npt) * 1000) / 0.002
time_smd = (float(args.smd) * 1000) / 0.002

def EnergyMinimization():
  # Create the system
  system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds, removeCMMotion=True)

  # Define the integrator
  integrator = LangevinIntegrator(float(args.temperature)*kelvin, 1.0/picosecond, 0.002*picoseconds)

  # Set up the simulation platform
  platform = Platform.getPlatformByName('CPU')
  properties = {"Threads":"6"}
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
  simulation.context.setPositions(inpcrd.positions)
  
  # Adding periodic condition
  if inpcrd.boxVectors is not None: 
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

  # Optimizing Protein Coordinates in the Simulation Box
  print("\t- Optimizing receptor coordinates along the x-axis..")
  protein_index = [atom.index for atom in pdb.topology.atoms() if atom.residue.chain.index in list(range(0,10))] 
  positions = simulation.context.getState(getPositions=True).getPositions()
  positions = np.array(positions / angstroms)
  positions[protein_index] -= [37, 0, 0]
  positions *= angstroms
  simulation.context.setPositions(positions)

  # Minimize the system's energy
  print("\t- Energy minimization started at {}".format(time.ctime()))
  simulation.minimizeEnergy()
  print("\t- Energy minimization finished at {}\n".format(time.ctime()))
  
  all_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
  PDBFile.writeFile(prmtop.topology, all_positions, open("em.pdb", 'w'))

def Equilibration():
  # Creating system and defining integrator
  system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)
  system.addForce(MonteCarloBarostat(1.0 * bar, float(args.temperature) * kelvin, 25))
  integrator = LangevinIntegrator(float(args.temperature) * kelvin, 1.0/picosecond, 0.002*picoseconds)
  pdb = PDBFile("em.pdb")
  
  # Applying harmonic position restraint on Ca atoms
  print("\t- Applying position restrain on protein Ca atoms")
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
  properties = {"Precision":"mixed",'DeviceIndex':str(args.gpuID)}
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
  simulation.context.setPositions(pdb.getPositions())
  simulation.context.setVelocitiesToTemperature(float(args.temperature) * kelvin)
  
  # Adding periodic condition
  if inpcrd.boxVectors is not None: 
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  
  # Running equilibration
  simulation.reporters.append(StateDataReporter("equilibration.log", 500, step=True, kineticEnergy=True, potentialEnergy=True,
    totalEnergy=True, temperature=True, volume=True))
  print("\t- Equilibration for {}ns started at {}".format(args.npt,time.ctime()))
  simulation.step(time_npt)
  print("\t- Equilibration finished at {}\n".format(time.ctime()))

  # Save rst7 and pdb file
  state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True, enforcePeriodicBox=True)
  with open("equilibration.rst7", 'w') as f:
    f.write(XmlSerializer.serialize(state))
    
  # Save PDB file
  all_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
  PDBFile.writeFile(prmtop.topology, all_positions, open("equilibration.pdb", 'w'))
  
def SMD():
  pdb = PDBFile("equilibration.pdb")
  system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds, removeCMMotion=True)
  integrator = LangevinIntegrator(float(args.temperature)*kelvin, 1.0/picosecond, 0.002*picoseconds)
  system.addForce(MonteCarloBarostat(1.0 * bar, float(args.temperature) * kelvin, 25))
  
  # Steered Molecular Dynamics (SMD) Setup
  plumed_script = """
  UNITS ENERGY=kcal/mol
  R1: COM ATOMS=5,15,35,49,65,79,103,109,131,143,162,183,199,215,230,251,258,269,283,300,314,333,348,359,381,409,415,431,446,468,483,502,514,533,543,553,572,591,607,628,652,667,684,699,711,733,747,766,785,802,822,838,855,862,877,892,904,923,945,961,978,995,1006,1017,1038,1062,1079,1103,1113,1137,1156,1175,1197,1209,1226,1245,1256,1275,1282,1296,1306,1316,1335,1352,1371,1385,1397,1413,1435,1454,1471,1483,1493,1500,1516,1537,1561,1572,1589,1608,1619,1640,1647,1654,1664,1676,1697,1719,1743,1762,1776,1792,1814,1830,1844,1862
  R2: COM ATOMS=1871,1882,1892,1902,1917,1927,1937,1961,1983,1993,2003,2022,2037,2059,2073,2083,2093,2104,2123,2138,2152,2162,2181,2191,2213,2223,2245,2260,2275,2282,2294,2309,2323,2339,2356,2378,2392,2411,2425,2435,2450,2469,2484,2506,2522,2546,2556,2571,2590,2605
  d: DISTANCE ATOMS=R1,R2 COMPONENTS
  x_dist: COMBINE ARG=d.x COEFFICIENTS=1 PERIODIC=NO
  MOVINGRESTRAINT ...
    LABEL=x_restraint
    ARG=x_dist
    STEP0=0 AT0=0.3 KAPPA0=0.0
    STEP1=250000 AT1=6.5 KAPPA1=250

  ... MOVINGRESTRAINT
  PRINT STRIDE=500 ARG=* FILE=COLVAR
  """
  force = PlumedForce(plumed_script)
  system.addForce(force)
  
  # Adding harmonic force restraint on CH3 region Ca atoms to avoid drifting during SMD
  restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
  restraint.addGlobalParameter('k', 50.0 * kilocalories_per_mole / nanometers**2)
  restraint.addPerParticleParameter('x0')
  restraint.addPerParticleParameter('y0')
  restraint.addPerParticleParameter('z0')
  
  for atom in pdb.topology.atoms():
      if 1 <= atom.residue.index <= 116 and atom.name == "CA":
          restraint.addParticle(atom.index, pdb.positions[atom.index])
  
  system.addForce(restraint)
  
  platform = Platform.getPlatformByName('CUDA')
  properties = {"Precision": "mixed", 'DeviceIndex': str(args.gpuID)}
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
  
  if inpcrd.boxVectors is not None:  # This will check if the system has periodic boundary condition
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  
  # Reading system state and coordinates
  with open("equilibration.rst7", 'r') as f:
      state = XmlSerializer.deserialize(f.read())
  simulation.context.setState(state)
  simulation.context.setPositions(state.getPositions())
  
  # Run SMD simulation
  print('\t- Running SMD...')
  simulation.reporters.append(DCDReporter('SMD.dcd', 500))
  simulation.reporters.append(PDBReporter('frames.pdb', 500,enforcePeriodicBox=True))
  simulation.reporters.append(StateDataReporter(stdout, 500, step=True, potentialEnergy=True, temperature=True))
  simulation.step(250000)
  # Save each frame
  
  # Saving simulation state
  state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True, enforcePeriodicBox=True)
  with open("smd.rst7", 'w') as f:
      f.write(XmlSerializer.serialize(state))

  # Saving final coordinates
  PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(), open("smd_output.pdb", 'w'))

if __name__ == "__main__":
  print("- Parsing coordinate file\n")
  inpcrd = AmberInpcrdFile("protein.inpcrd")
  print("- Parsing topology file\n")
  prmtop = AmberPrmtopFile("protein.prmtop")
  print("- Parsing pdb file\n")
  pdb = PDBFile("protein_WAT.pdb")
  
  print("- Running energy minimization\n")
  #EnergyMinimization()
  
  print("- Running equilibration\n")
  Equilibration()
  
  print("- Running SMD\n")
  SMD()
  
  print('SMD Simulation is Completed!')
