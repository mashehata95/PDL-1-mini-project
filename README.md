# Protein Binder Design for PDL-1

## Introduction
Programmed death-ligand 1 (PD-L1) is an immune checkpoint protein that plays a crucial role in downregulating immune responses, allowing cancer cells to evade immune surveillance. Inhibiting PD-L1 has been a successful strategy in cancer immunotherapy, with monoclonal antibodies like Atezolizumab and Durvalumab demonstrating significant clinical benefits. However, the design of novel protein-based binders for PD-L1 provides an alternative therapeutic approach with potential advantages in specificity, stability, and manufacturability.

This mini-project focuses on leveraging computational protein design techniques to develop a novel protein binder targeting PD-L1. By integrating molecular dynamics (MD) simulations, protein interaction fingerprint (PIF) analysis, and deep learning-based protein design using RFdiffusion, this project aims to generate and validate high-affinity PD-L1 binders.

## Project Outline
This project follows a structured computational pipeline:

1. **Acquisition of Experimental Structures**
   - Obtain the experimentally resolved structures of Atezolizumab-PD-L1 (PDB ID 5X8L) and Durvalumab-PD-L1 (PDB ID 5X8M) from the Protein Data Bank (PDB).

2. **Structure Processing & Preparation for MD Simulations**
   - Preprocess the structures by removing irrelevant components (e.g., solvent, ions) and optimizing for molecular dynamics (MD) simulations.
   - Generate topology and parameter files.

3. **100 Nanoseconds MD Simulations**
   - Perform classical MD simulations to study the conformational flexibility and interaction dynamics of the PD-L1-antibody complexes using OpenMM.

4. **Protein Interaction Fingerprint (PIF) Analysis**
   - Identify key residues at the PD-L1 interface involved in critical interactions with Atezolizumab and Durvalumab.
   - Define target hotspots for designing novel protein binders.

5. **RFdiffusion-based Protein Binder Design**
   - Utilize RFdiffusion to design novel protein binders specifically targeting the identified hotspot residues. The design was run on the following colab: https://colab.research.google.com/github/sokrypton/ColabDesign/blob/v1.1.1/rf/examples/diffusion.ipynb

6. **Selection of Top Binders Based on RMSD Profile**
   - Filter generated binders based on their structural deviation and stability compared to reference structures.

7. **MD Simulation of PD-L1 and Selected Binders**
   - Conduct additional MD simulations to confirm the stability of the PD-L1-binder complexes.

8. **Extraction of Relaxed Conformation**
   - Identify the most stable binding conformations from MD trajectories.

9. **Free Energy Evaluation Using Umbrella Sampling**
   - Compute the binding free energy (ΔG) of the designed binders to assess their binding affinity and stability.

## Repository Contents
This repository contains:

- **Notebooks/**
  - `structure_processing.ipynb`: Prepares experimental structures for MD simulations.
  - `fingerprint_analysis.ipynb`: Identifies PD-L1 interaction hotspots using Protein Interaction Fingerprinting.

- **Results/**
  - `presentation.pdf`: Highlights key findings from the project, including interaction analysis, designed binders, and stability evaluation.


For detailed instructions, refer to the documentation in the respective notebooks.

## Contact
For questions or collaborations, feel free to reach out to **mashehata95@gmail.com**.


