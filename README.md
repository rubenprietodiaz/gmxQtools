# gmxQTools
This repository provides tools and scripts for working with molecular dynamics (MD) simulations using GROMACS (e.g., PyMemDyn), as well as for performing simulations with Q (e.g., QligFEP) on the same molecular system. It includes tools for system preparation, simulation execution, results analysis, and more.

All the scripts and the pipeline are optimized for batch processing of MD simulations, handling numerous ligands, folders, and proteins.

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

### 1.2 Protein and Ligand Preparation

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
   
    ```bash
        setup_pym [-C CLUSTER] [-p PROTEIN]
                [-l LIGAND] [-w WATERS]
                [-i IONS] [--full_relax FULL_RELAX]
                [--res RESTRAINT]
        
        -h, --help
                    show help message

        -C CLUSTER
                    Choose your cluster over the list.
                    You can add more by modifying the code.
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

        --full_relax [True/False]
                    Toggle for performing full MD relaxation.
                    If set to false, the run will finish after
                    the initial relaxation (default = True)
                    See section 2 to choose.
    ```

   
This creates a folder for each ligand, executes ligpargen for ligand parameters and generate scripts for pymemdyn execution (pymemdyn.sh inside ligand folder and submit.sh for        submitting to a cluster SLURM queue.

### 2. Execute pymemdyn
Execute submit.sh file in your cluster.
Check the original repository [PyMemDyn](https://github.com/GPCR-ModSim/pymemdyn) for requirements, installation and tutorials. The necessary arguments were included in submit.sh script in step 1.2.6, but you can modify pymemdyn.sh inside each folder with your preferences, or the generation of this script in `pym_setup.py`.

From this point onward, the protocol will vary depending on whether you intend to perform MD simulations, FEP calculations, or both. Choose the appropriate protocol based on your objectives.

#### 2.1 Free Energy Perturbation (FEP) with QligFEP
If you do not have any intention of running MD simulations, it is recommended to choose `--full_relax False` when executing `setup_pym.py`. After that, choose option -n 1 to prepare your system for FEP when running `setup_md.py`. (TO BE INCLUDED)

#### 2.2 Molecular Dynamics (MD) or both
If you want to do MD simulations or both MD and FEP, choose choose `--full_relax True` (or leave it blank) when executing `setup_pym.py`. After that, choose option -n 2 to prepare your system for FEP when running `setup_md.py`. (TO BE INCLUDED)

 ```bash
      setup_md  [-C CLUSTER] [-n OPTION]
                [-d DIR] [-t TIME]
      
      Optional arguments (for executing pymemdyn after this preparation:

      -h, --help
                    show help message

      -C CLUSTER
                    Choose your cluster over the list.
                    You can add more by modifying the code.
      
      -n OPTION
                    Choose between prepare the system for:
                    -n 1 for FEP
                    -n 2 for FEP and/or MD
      
      -d DIR
                    Directory for input md files (default:
                    md_input_files)
      -t TIME (ns)
                    Time for MD simulation (in nanoseconds)
 ```


