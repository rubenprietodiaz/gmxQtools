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
rpmol {file_to_convert}
```

The source code of rpmol is also included in '$FOLDER'.
