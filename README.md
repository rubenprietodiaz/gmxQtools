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
All the required versions are included in ´requirements.txt´.

## Setup Scripts Overview

### 1. `setup_pym.py`
This script prepares the system for PyMemDyn simulations. It creates complexes between ligand and receptor, aligns the protein, generates ligand parameters using Ligpargen, and renames the files properly for PyMemDyn.

#### Usage

```bash
setup_pym [-C CLUSTER] [-p PROTEIN]
          [-l LIGAND] [-w WATERS]
          [-i IONS] 
          [--res RESTRAINT]
```

- **-h, --help**: Show help message
- **-C CLUSTER**: Choose your cluster from the list (modify code to add more clusters)
- **-p PROTEIN**: PDB file of your protein (default = protein.pdb)
- **-l LIGAND**: Ligand identifiers in pdb (default = UNK, from Maestro)
- **-w WATERS**: Water identifiers of crystallized water molecules in the PDB file
- **-i IONS**: Ion identifiers of crystallized ions in the PDB file
- **--res RESTRAINT**: Position restraints during MD production run. Options: bw (Ballesteros-Weinstein Restrained Relaxation - default), ca (C-Alpha Restrained Relaxation)
- **--fep**: If selected, the run will finish after the initial relaxation. Choose this option if you know you only want to perform FEP simulations.

After execution, this script creates a folder for each ligand, executes Ligpargen for ligand parameters, and generates scripts for PyMemDyn execution (`pymemdyn.sh` inside the ligand folder and `submit_pym.sh` for submitting to a cluster SLURM queue).

### 2. `setup_md.py`
This script prepares the files for running an MD simulation. It should be executed inside your directory containing subdirectories for each ligand after running PyMemDyn.

#### Usage

```bash
setup_md [-C CLUSTER] [-t TIME]
         [-rt RUNTIME]
```

- **-h, --help**: Show help message
- **-C CLUSTER**: Choose your cluster from the list (modify code to add more clusters)
- **-t TIME (ns)**: Time for MD simulation (in nanoseconds)
- **-rt RUNTIME (hours)**: Limit of time for simulation (in hours)

After execution, the script generates `md_input_files`. To submit the MD job, enter the directory and execute `sh submit_md.sh`.

### 3. `setup_fep.py`
This script prepares files for FEP simulations. It should be executed inside your directory containing subdirectories for each ligand after running `setup_pym.py`.

#### Usage

```bash
setup_fep [-d DIR] [-nc]
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

## Analysis Scripts Overview
Under construction.

## Example GPCR Workflow
<p align="center">
  <img src="/manual/Protocol.jpg" alt="Protocol" />
</p>
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

    Please, note that PyModSim execution is included in setup_pym protocol, through: `pymodsim -n 3 -p [PDB]`, which means **no corrections** in the structure of the PDB file. Check [Pymodsim](https://github.com/GPCR-ModSim/pymodsim) for more info. If you want to perform some corrections, after preparing the protein in Maestro, run PyModSim manually and continue the protocol with `old_setup_pym.py` located in developing directory.

3. **Submit PyMemDyn Jobs**:
    - Execute `sh submit_pym.sh` to submit PyMemDyn jobs to your cluster.

4. **Prepare MD Simulation**:
    - Run `setup_md.py` to generate MD input files.
    - This will create a directory for each ligand, containing the necessary input files for running the MD simulations.
    - Submit the MD job by executing `sh submit_md.sh`.
    - Analyze data with analysis scripts provided.

5. **Prepare FEP Simulation** (if needed):
    - Run `setup_fep.py` to prepare FEP files.
    - This will create the directory fep_preparation_files that you can use for modelling novel ligands and running QligFEP protocol.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [PyMemDyn](https://github.com/GPCR-ModSim/pymemdyn) for providing the tools for membrane molecular dynamics calculations.
- [PyModSim](https://github.com/GPCR-ModSim/pymodsim) for providing the tools for membrane alignment.
- [Ligpargen](https://github.com/Isra3l/ligpargen) for ligand parameter generation.

For any issues or contributions, please open a ticket or submit a pull request on the project's GitHub page.
