# gmxQTools
This repository provides tools and scripts for working with molecular dynamics (MD) simulations using GROMACS (e.g., PyMemDyn), as well as for performing simulations with Q (e.g., QligFEP) on the same molecular system. It includes tools for system preparation, simulation execution, results analysis, and more.

All the scripts and the pipeline are optimized for batch processing of MD simulations, handling numerous ligands, folders, and proteins.


## Installation
Clone the repository to your machine:

```bash
git clone https://github.com/your_username/gmxQTools.git
```

Ensure that you have the necessary dependencies installed, including GROMACS, Q, QligFEP, PyMemDyn, PyModSim, Ligpargen, Modeller, and Python packages specified in the `requirements.txt` file.

## GPCR Work Pipeline

This pipeline outlines the workflow for working with G Protein-Coupled Receptors (GPCRs). All steps described below have been tested using Schrödinger software (Maestro) for file preparation (protprep, ligprep), and docking.

### 1. File preparation

#### 1.1 Excel to SDF Conversion

Quickly convert an Excel file containing IDs and SMILES to SDF format using `rpmol`. This can be done easily on your computer with a simple installation via pip:

```bash
pip install rpmol
```

After installation, you can convert xlsx files containing 'ID' and 'SMILES' columns to sdf. Also allows to convert sdf files to xlsx.

```bash
rpmol [file_to_convert]
```

The source code of rpmol is also included in this repository.

#### 1.2 Protein and Ligand Preparation

In Maestro, protein and ligand preparation typically involve several steps to ensure proper structure and compatibility for molecular docking studies. Here's a brief overview:

1. **Protein Preparation**:
    - Import the protein structure (usually in PDB format).
    - Execute Protein Preparation tool in Maestro.

2. **Ligand Preparation**:
    - Import the ligand structures (sdf file generated before).
    - Execute Ligand Preparation tool in Maestro.

3. **Pymodsim centering**: As the GPCR should be embedded in a membrane, you first need to align the protein with it.
    - Export prepared protein as pdb file.
    - Follow the protocol of [PyModSim](https://github.com/GPCR-ModSim/pymodsim). For only alignment:

      ```bash
      pymodsim -n 3 -p [PDB]
      ```
      
4. **Grid generation and ligand docking**: Import pymodsim aligned protein and create a grid to your specific requirements. Then proceed with ligand docking.
   
5. **Export files**:
    - Export selected poses as pdb files.
    - Export GPCR as pdb file. We recommend to name it protein.pdb (default for -p in next step).

6. **Prepare system for PyMemDyn**: at this point you have multiple ligands (and/or diferent poses for same ligand) as pdb files with residue name 'UNK' and the protein.

Execute 'setup_pym.py' in the directory containing all the files to create complexes between ligand and receptor, generate parameters of the ligands using [Ligpargen](https://github.com/Isra3l/ligpargen), and rename the files properly for PyMemDyn. 
   
        setup_pym [-C CLUSTER] [-p PROTEIN]
                  [-l LIGAND] [-w WATERS]
                  [-i IONS] 
                  [--res RESTRAINT]
        
        -h, --help
                    show help message

        -C CLUSTER
                    Choose your cluster over the list.
                    You can add more by modifying the code
        -p PROTEIN 
                    PDB file of your protein
                    (default = protein.pdb)
        -l LIGAND
                    Ligand identifiers of ligand in pdb.
                    (default = UNK, from Maestro)

    Optional arguments for executing pymemdyn after this preparation:

        
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

        --fep 
                    If selected, the run will finish after
                    the initial relaxation.

   
This creates a folder for each ligand, executes ligpargen for ligand parameters and generate scripts for pymemdyn execution (pymemdyn.sh inside ligand folder and submit_pym.sh for submitting to a cluster SLURM queue.

### 2. PyMemDyn execution
PyMemDyn is a standalone python package to setup membrane molecular dynamics calculations using the GROMACS set of programs. In the previous steps, you have prepared your system for PyMemDyn. For running it, execute submit_pym.sh file in your cluster.

```bash
sh submit_pym.sh
```

Check the original repository [PyMemDyn](https://github.com/GPCR-ModSim/pymemdyn) for requirements, installation and tutorials. The necessary arguments were included in submit_pym.sh script in step 1.2.6, but you can modify pymemdyn.sh inside each folder with your preferences, or the generation of this script in `pym_setup.py`.

### 3. Preparation of MD
For preparing the files for running an MD simulation, execute `setup_md.py` inside your directory (containing subdirectories for each ligand) after running pymemdyn.

```bash
        setup_md  [-C CLUSTER] [-t TIME]
                  [-rt RUNTIME]
        
        -h, --help
                    show help message

        -C CLUSTER
                    Choose your cluster over the list.
                    You can add more by modifying the code
        
        -t TIME (ns)
                    Time for MD simulation (in nanoseconds)

        -rt RUNTIME (hours)
                    Limit of time for simulation (in hours)
```

After that, you will get md_input_files. For submitting MD, to enter and execute `sh submit_md.sh`.

### 4. Preparation of FEP
For preparing files for FEP, execute `setup_fep.py` inside your directory (containing subdirectories for each ligand) after running pymemdyn.

```bash
        setup_fep [-d DIR] [-nc]
        
        -h, --help
                    show help message

        -d DIR
                    Choose your directory for creating FEP
                    files.

        -nc, --noclean
                    Keep temporary files and log in the output.
```

After that, you will get fep_preparation_files: 
    
    - **system.pdb**: Trimmed system containing membrane, water, and protein for FEP, with solvent removal and changes in nomenclature.
    - **ligand.pdb**: Template for model different ligands for FEP.
    
If using --noclean argument:

    - **complex.pdb**: PDB file of the unprepared system.
    - **water.pdb**: Extracted waters unprepared.
    - **water_fixed.pdb**: Renumbered waters.
    - **pymol.log**: Log of pymol transformations.
    - **system.pdb.bak**: Backup of system.pdb before renaming residues for Q.
    - **trjconv.log**: Output of gmx trjconv script to convert *.gro to *.pdb.



