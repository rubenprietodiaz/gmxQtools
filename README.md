# gmxQTools
This repository provides tools and scripts for working with molecular dynamics (MD) simulations using GROMACS (e.g., PyMemDyn), as well as for performing simulations with Q (e.g., QligFEP) on the same molecular system. It includes tools for system preparation, simulation execution, results analysis, and more.

## GPCR Work Pipeline

This pipeline outlines the workflow for working with G Protein-Coupled Receptors (GPCRs). All steps described below have been tested using Schrödinger software (Maestro) for file preparation (protprep, ligprep), and docking.

### 1. Protein and Ligand Preparation

#### 1.1 Excel to SDF Conversion

Quickly convert an Excel file containing IDs and SMILES to SDF format using `rpmol`. This can be done easily on your computer with a simple installation via pip:

```bash
pip install rpmol
```
After installation, you can convert xlsx files containing 'ID' and 'SMILES' columns to sdf. Also allows to convert sdf files to xlsx.
```bash
rpmol [file_to_convert]
```

The source code of rpmol is also included in '$FOLDER'.

### 1.2 Protein and Ligand Preparation

In Maestro, protein and ligand preparation typically involve several steps to ensure proper structure and compatibility for molecular docking studies. Here's a brief overview:

1. **Protein Preparation**:
    - Import the protein structure (usually in PDB format).
    - Execute Protein Preparation tool in Maestro.

2. **Ligand Preparation**:
    - Import the ligand structures (sdf file generated before)
    - Execute Ligand Preparation tool in Maestro.

3. **Pymodsim centering**: As the GPCR will be embedded in a membrane, you first need to align the protein with it.
    - Export prepared protein as pdb file.
    - Follow the protocol of [PyModSim](https://github.com/GPCR-ModSim/pymodsim). For only alignment:
      ```bash
      pymodsim -n 3 -p [PDB]
      ```
      
4. **Grid generation and ligand docking**: Create a grid to your specific requirements, and then proceed with ligand docking.
   
5. **Export files**:
    - Export selected poses as pdb files.
    - Export GPCR as pdb file, save as protein.pdb.

6. **Prepare system for pymemdyn**: at this point you have multiple ligands (and/or diferent poses for same ligand) as pdb files with residue name 'UNK' and the protein.
   Execute 'setup_pym.py' in the directory containing all the files to create complexes between ligand and receptor, generate parameters of the ligands using [Ligpargen][https://github.com/Isra3l/ligpargen], and rename the files properly for pymemdyn.
      ```bash
      setup_pym [-C CLUSTER] [-l LIGAND]
                [-w WATERS] [-i IONS]
                [--full_relax FULL_RELAX]
                [--res RESTRAINT]
      Optional arguments (for executing pymemdyn after this preparation:
      -h, --help
                    show help message
      -l LIGAND
                    Ligand identifiers of ligands
                    present within the PDB file.
                    If multiple ligands are present,
                    give a comma-delimited list.
      
      -w WATERS
                    Water identifiers of crystalized
                    water molecules present within
                    the PDB file.
      -i IONS
                    Ion identifiers of crystalized water
                    molecules present within the PDB file.
      
      --res [RESTRAINT]
                    Position restraints during MD production
                    run. Options: bw (Ballesteros-Weinstein
                    Restrained Relaxation - default),
                    ca (C-Alpha Restrained Relaxation)

      --full_relax [True/False]
                    Toggle for performing full MD relaxation.
                    If set to false, the run will finish after
                    the initial relaxation. (default = True)
      ```

Edit the code for adding your own clusters.
    
This creates a folder for each ligand, executes ligpargen for ligand parameters and generate scripts for pymemdyn execution (pymemdyn.sh inside ligand folder and submit.sh for           submitting to a cluster SLURM queue.

7. **Run pymemdyn**: execute submit.sh file in your cluster.




These steps ensure that both the protein and ligands are in a suitable state for subsequent docking studies. For detailed protocols and specific tools within Maestro, refer to the Schrödinger documentation or tutorials.
