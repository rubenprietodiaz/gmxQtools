# MD and FEP setup

This repository contains scripts for setting up molecular dynamics (MD) simulations using GROMACS and free energy perturbation (FEP) simulations using Q for ligand-GPCR complexes. The provided scripts simplify the preparation process for these simulations, making it easier to configure and execute them using the pymemdyn package. Additionally, these scripts allow you to prepare multiple ligands simultaneously.

## Prerequisites

Before using these scripts, ensure you have the following software installed:
- Python (3.7 or higher): highly recommended to create a conda environment named py37.
- `PyModSim`: install via [GitHub](https://github.com/GPCR-ModSim/pymodsim).
- `PyMemDyn`: install via [GitHub](https://github.com/GPCR-ModSim/pymemdyn).
- `Ligpargen`: install via [GitHub](https://github.com/Isra3l/ligpargen).
- Schrödinger Maestro: Required for protein and ligand preparation. Note: These scripts have only been tested with the Schrödinger Suite. Alternatively, you may use similar software if compatible.
- GROMACS: required for MD simulations.
- PyMol: for setup FEP files via `conda install -c conda-forge pymol-open-source`. 

For analysis:
- mdtraj
- pandas
- numpy

All the required modules are included in requirements.txt. You may also use the environment.yml file in this repository. Ensure the installation of LigParGen, PyModSim, and PyMemDyn by following the instructions located in each respective repository.

```bash
conda env create -f environment.yml
conda activate py37
```

## 1. Setup Scripts Overview

### A. `setup_pym.py`
This script prepares the system for PyMemDyn simulations. It creates complexes between ligand and receptor, aligns the protein, generates ligand parameters using Ligpargen, and renames the files properly for PyMemDyn.

#### Usage

```bash
setup_pym.py [-C CLUSTER] [-p PROTEIN]
             [--noclean] [--noalign]
             [-cb CHARGE] [-l LIGAND]
             [-w WATERS] [-fep] 
             [-i IONS] [-r RESTRAINT]             
```

- **-h, --help**: Show help message
- **-C CLUSTER**: Choose your cluster from the list (modify code to add more clusters)
- **-p PROTEIN**: PDB file of your protein (default = protein.pdb)
- **-l LIGAND**: Ligand identifiers in pdb (default = UNK, from Maestro)
- **-w WATERS**: Water identifiers of crystallized water molecules in the PDB file
- **-i IONS**: Ion identifiers of crystallized ions in the PDB file
- **-cb CHARGE**: Ligand charge (for ligpargen execution, default = 0)
- **-r RESTRAINT**: Position restraints during MD production run. Options: bw (Ballesteros-Weinstein Restrained Relaxation - default), ca (C-Alpha Restrained Relaxation)
- **--fep**: If selected, the run will finish after the initial relaxation. Choose this option if you know you only want to perform FEP simulations.
- **--noclean**: If selected, the directory will not be clean after processing.
- **--noalign**: Use this option if the complex is already aligned with PyModSim.

After execution, this script creates a folder for each ligand, executes Ligpargen for ligand parameters, and generates scripts for PyMemDyn execution (`pymemdyn.sh` inside the ligand folder and `submit_pym.sh` for submitting to a cluster SLURM queue).

**Before execution**:
```bash
my_project
├── protein.pdb
├── lig1.pdb
├── lig2.pdb
├── lig3.pdb
├── lig(n).pdb
└── setup_pym.py <- execute this to prepare files
```

**After setup_pym.py execution**:
```bash
my_project
├── inputFiles
│   ├── protein.pdb
│   ├── lig1.pdb
│   ├── lig2.pdb
│   ├── lig3.pdb
│   └── lig(n).pdb
├── lig1
│   ├── complex.pdb <- Aligned complex between lig1.pdb and protein.pdb
│   ├── files for ligand parameters (*.itp, *.rtf, *.cms...)
│   ├── ligpargen.log <-ligpargen execution log
│   └── pymemdyn.sh <- SLURM execution for PyMemDyn
├── lig2
├── lig3
├── lig(n)
└── submit_pym.sh  <- Calls pymemdyn.sh inside each directory

```

**After running PyMemDyn**:
```bash
my_project
├── inputFiles
│   ├── protein.pdb
│   ├── lig1.pdb
│   ├── lig2.pdb
│   ├── lig3.pdb
│   └── lig(n).pdb
├── lig1
│   ├── files & directories
│   ├── eqProd/
│   │   ├── confout200.gro <- Starting point for FEP
│   │   ├── confout.gro
│       └── other files
│   └── finalOutput
│       ├── logs/
│       ├── reports/
│       ├── prod.mdp <- Template for MD
│       ├── confout.gro <- Starting point for MD
│       └── module2
├── lig2
├── lig3
├── lig(n)
└── submit_pym.sh  

```
### B. `setup_md.py`
This script prepares the files for running an MD simulation. It should be executed inside your directory containing subdirectories for each ligand after running PyMemDyn.

#### Usage

```bash
setup_md.py [-C CLUSTER] [-t TIME]
            [-rt RUNTIME] [-n NUMREPLICAS]
```

- **-h, --help**: Show help message
- **-C CLUSTER**: Choose your cluster from the list (modify code to add more clusters)
- **-t TIME (ns)**: Time for MD simulation (in nanoseconds)
- **-rt RUNTIME (hours)**: Limit of time for simulation (in hours)
- **-n NUM_REPLICAS**: Number of replicas for the MD simulations (default: 1).

After execution, the script generates `md_input_files`. To submit the MD job, enter the directory and execute `sh submit_md.sh`.

**After running setup_md.py**:
```bash
my_project
├── inputFiles
│   ├── protein.pdb
│   ├── lig1.pdb
│   ├── lig2.pdb
│   ├── lig3.pdb
│   └── lig(n).pdb
├── lig1
├── lig2
├── lig3
├── lig(n)
├── submit_pym.sh
├── setup_md.py
└── md_input_files <- New folder
│   ├── lig1_1 <- Replica 1 ligand 1
│   │   ├── LIG.itp         
│   │   ├── ffoplsaa_mod.itp
│   │   ├── index.ndx  
│   │   ├── posre.itp      
│   │   ├── posre_NA.itp  
│   │   ├── run_md.sh <- SLURM execution for MD
│   │   ├── topol.tpr
│   │   ├── LIG_backup.itp  
│   │   ├── ffoplsaabon_mod.itp  
│   │   ├── ions.itp   
│   │   ├── posre_HOH.itp
│   │   ├── prod.mdp <- Modified with setup_md.py
│   │   ├── spc.itp
│   │   ├── confout.gro
│   │   ├── ffoplsaanb_mod.itp
│   │   ├── popc.itp
│   │   ├── posre_LIG.itp
│   │   ├── protein.itp
│   │   ├── topol.top
│   ├── lig1_2 <- Replica 2 ligand 1
│   ├── lig1_n <- Replica n ligand 1
│   └── submit_md.sh <- Calls run_md.sh inside each directory

```

### C. `setup_fep.py`
This script prepares files for FEP simulations. It should be executed inside your directory containing subdirectories for each ligand after running `setup_pym.py`.

#### Usage

```bash
setup_fep.py [-d DIR] [-nc]
```

- **-h, --help**: Show help message
- **-d DIR**: Choose your directory for creating FEP files, taking this ligand as a model.
- **-nc, --noclean**: Keep temporary files and log in the output

After execution, the script generates `fep_preparation_files`:
- **system.pdb**: Trimmed system containing membrane, water, and protein for FEP, with solvent removal and changes in nomenclature
- **ligand.pdb**: Template for modeling different ligands for FEP

If using the `--noclean` argument:
- **complex.pdb**: PDB file of the unprepared system
- **water.pdb**: Extracted waters unprepared
- **water_fixed.pdb**: Renumbered waters
- **pymol.log**: Log of pymol transformations
- **system.pdb.bak**: Backup of system.pdb before renaming residues for Q
- **trjconv.log**: Output of gmx trjconv script to convert *.gro to *.pdb

**After running setup_md.py**:
```bash
my_project
├── inputFiles
│   ├── protein.pdb
│   ├── lig1.pdb
│   ├── lig2.pdb
│   ├── lig3.pdb
│   └── lig(n).pdb
├── lig1
├── lig2
├── lig3
├── lig(n)
├── submit_pym.sh
├── setup_md.py
└── fep_preparation_files <- New folder
│   ├── system.pdb <- Membrane, cofactors (ions), solvent and protein.
│   └── ligand.pdb <- Ligand for modelling
```


## 2. Analysis Scripts Overview
### A. `rmsd_lig.py`
This script calculates the RMSD of a ligand in trajectory files, grouping RMSD values from replicated assays into a single file. It processes directories, computes RMSD values, and saves the results in both .xvg (frame by frame) and .txt (mean in each replicate) formats. The script now supports SMARTS patterns and SMILES codes for analyzing specific parts of the ligand, and the option to calculate RMSD for atoms not matching the SMARTS pattern or MCS with SMILES.

#### Usage

```bash
rmsd_lig.py [-p PDB_FILE] [-t TRAJ_FILE]
            [-s SMARTS_PATTERN] [-S SMILES]
            [-i] [-l LIGAND_NAME]
            [-o OUTPUT_FILENAME] 
            [-f REFERENCE_FRAME]
```

- **-h, --help**: Show help message.
- **-p PDB_FILE**: Path to the PDB file (default: `finalOutput/start.pdb`).
- **-s SMARTS_PATTERN**: Deprecated. SMARTS pattern for specific part of the ligand. Provide a simplified pattern (no Hs, no aromaticity, no bond orders).
- **-S SMILES**: SMILES code. If provided, calculate maximum common substructure (MCS) with the SMILES. Optionally analyze a specific part of the ligand by providing the SMILES code for that part.
- **-i, --inverse**: Calculate RMSD for atoms not matching the SMARTS pattern or MCS with SMILES.
- **-l LIGAND_NAME**: Name of the ligand (default: `L01`).
- **-t TRAJ_FILE**: Path to the trajectory file (default: `finalOutput/traj_prod_pymol.xtc`).
- **-o OUTPUT_FILENAME_RMSD**: Output filename for RMSD results (default: `rmsd_stat.txt`).
- **-f REFERENCE_FRAME**: Frame to use as reference for RMSD calculation (default: 0).

#### Functionality

- **SMARTS Pattern**: Analyze specific parts of the ligand by providing a SMARTS pattern. Deprecated, and users are encouraged to use SMILES for better accuracy.
- **SMILES Code**: Calculate maximum common substructure (MCS) with the ligand using a provided SMILES code. This also allows for analysis of specific parts of the ligand.
- **Inverse Selection**: Calculate RMSD for atoms that do not match the SMARTS pattern or MCS with SMILES by using the `-i` flag.
- **Exclusion of Hydrogens**: Hydrogens are excluded from RMSD calculations to ensure consistency and relevance of results.
- **Combined RMSD Values**: Saves combined RMSD values for each group to a single .xvg file, facilitating easier analysis of replicated assays.

## Example GPCR setup Workflow
1. **Prepare Protein and Ligand**:
    - Use Schrödinger Maestro to prepare protein and ligand structures.
    - Export prepared protein as `protein.pdb` and ligands as PDB files with residue name 'UNK'. 
    - **Highly recommended**: It is highly recommended to use the default file type names throughout the protocol.
    - **Note**: Although Maestro allows you to change the residue name from UNK to anything you want, the `setup_pym` script is optimized for renaming from UNK.

2. **Execute `setup_pym.py`**:
    - Place all files (protein and ligands) in a directory.
    - Run `setup_pym.py` to prepare the system for PyMemDyn.
    - This creates a subdirectory for each ligand, containing:
        - The ligand-receptor complex file, aligned with the membrane (through PyModSim) for properly PyMemDyn execution.
        - Ligand parameters generated by Ligpargen.
        - A script (`pymemdyn.sh`) to execute PyMemDyn within each ligand's subdirectory.
        - A submission script (`submit_pym.sh`) to submit the PyMemDyn jobs to a cluster SLURM queue (by calling to `pymemdyn.sh`.).

    Please, note that PyModSim execution is included in setup_pym protocol, through: `pymodsim -n 3 -p [PDB]`, which means **no corrections** in the structure of the PDB file. Check [Pymodsim](https://github.com/GPCR-ModSim/pymodsim) for more info. If you want to perform some corrections, after preparing the protein in Maestro, run PyModSim manually and continue the protocol with `--noalign` argument.

3. **Submit PyMemDyn Jobs**:
    - Execute `sh submit_pym.sh` to submit PyMemDyn jobs to your cluster.

4. **Prepare MD Simulation**:
    - Run `setup_md.py` to generate MD input files.
    - This will create a directory for each ligand, containing the necessary input files for running the MD simulations.    
    - Submit the MD job by executing `sh submit_md.sh`.
    - Analyze data with analysis scripts provided *(Section 2. Analysis Scripts Overview)*.

5. **Prepare FEP Simulation** (if needed):
    - Run `setup_fep.py` to prepare FEP files.
    - This will create the directory fep_preparation_files that you can use for modelling novel ligands and running [QligFEP](https://github.com/qusers/qligfep) protocol.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [PyMemDyn](https://github.com/GPCR-ModSim/pymemdyn) for providing the tools for membrane molecular dynamics calculations.
- [PyModSim](https://github.com/GPCR-ModSim/pymodsim) for providing the tools for membrane alignment.
- [Ligpargen](https://github.com/Isra3l/ligpargen) for ligand parameter generation.
- @JarrettSJohnson for [error](https://github.com/schrodinger/pymol-open-source/issues/304) solution with Pymol.

For any issues or contributions, please open a ticket or submit a pull request on the project's GitHub page.
