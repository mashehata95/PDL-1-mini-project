# Umbrella Sampling Protocol

## Requirements
- OpenMM
- OpenMM-Plumed

## Overview
This protocol describes the steps to perform Umbrella Sampling to study the binding mechanism of a binder protein to PDL1.

## 1. Preprocessing
- The relaxed conformation of the complex was extracted as a function of its RMSD.
- The structure was rotated to orient the binder protein in the positive X-axis direction.
- The topology, coordinates, and solvated PDB files were generated using `tleap`.
- The system was solvated in a water box with dimensions of **60 x 20 x 20** Å to ensure sufficient space for the steering process.

## 2. Steered Molecular Dynamics (SMD) Simulation
- A steered MD simulation was performed to separate the binder from the PDL1 protein.
- The center of mass (COM) of two groups of atoms (R1 and R2) was calculated.
- The X-axis component of the COM distance was used as the reaction coordinate.
- A moving harmonic restraint was applied, starting from an initial separation of **0.3 Å** and increasing to **6.5 Å** over **250,000 simulation steps**.
- The force constant was ramped from **0.0 to 250 kcal/mol/Å²**.

## 3. Extraction of Umbrella Windows
- The steered trajectory was processed using the `SMD_processing.py` script.
- Frames with equidistant COM spacing of **1 Å** were extracted, generating **69 umbrella windows**.

## 4. Umbrella Sampling Simulations
Each window underwent the following steps:
- **Energy Minimization**: To remove steric clashes and optimize the structure.
- **Heating**: Conducted in the NVT ensemble at **310 K**, with strong harmonic restraints on Cα atoms.
- **Equilibration**: Performed in the NPT ensemble at **310 K**, with strong harmonic restraints on Cα atoms.
- **Production**:
  - The COM distance between R1 and R2 along the X-axis was used as the reaction coordinate.
  - A harmonic restraint of **10 kcal/mol/Å²** was applied at the COM distance of each window.
  - Simulation output was stored every **500 steps**.

## 5. Free Energy Analysis
- The Weighted Histogram Analysis Method (WHAM) was used to construct the **Potential of Mean Force (PMF) profile**.
- The **binding free energy (ΔG)** was computed from the PMF profile.

This protocol ensures a robust and systematic approach to studying the molecular interactions between the binder protein and PDL1 using umbrella sampling.
