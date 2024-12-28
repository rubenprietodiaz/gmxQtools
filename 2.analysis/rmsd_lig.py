import os
import argparse
import mdtraj as md  # type: ignore
import pandas as pd  # type: ignore
import numpy as np  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
from rdkit.Chem import rdFMCS  # type: ignore
from collections import defaultdict

# Argument parser setup and input parameters
parser = argparse.ArgumentParser(description="Calculate RMSD of a ligand in trajectory files.")
parser.add_argument('-p', '--pdb_file', type=str, default='finalOutput/start.pdb', help='Path to the PDB file. Default is finalOutput/start.pdb')
parser.add_argument('-s', '--smiles', type=str, help='SMILES code. If provided, calculate maximum common substructure (MCS) with the SMILES. If you want to analyze a specific part of the ligand, you can provide the SMILES code for that part.')
parser.add_argument('-i', '--inverse', action='store_true', help='Calculate RMSD for atoms not matching the MCS with SMILES')
parser.add_argument('-l', '--ligand_name', type=str, default='L01', help='Name of the ligand. Default is L01')
parser.add_argument('-t', '--traj_file', type=str, default='finalOutput/traj_prod_pymol.xtc', help='Path to the trajectory file. Default is finalOutput/traj_prod_pymol.xtc')
parser.add_argument('-o', '--output_filename_rmsd', type=str, default='rmsd_stat.txt', help='Output filename for RMSD results. Default is rmsd_stat.txt')
parser.add_argument('-f', '--reference_frame', type=int, default=0, help='Frame to use as reference for RMSD calculation. Default is 0')

args = parser.parse_args()

# Assign parsed arguments to variables
pdb_file = args.pdb_file
smiles = args.smiles
inverse = args.inverse
ligand_name = args.ligand_name
traj_file = args.traj_file
output_filename_rmsd = args.output_filename_rmsd
reference_frame = args.reference_frame

# Check for conflicting arguments
if inverse and not smiles:
    print("Inverse flag provided without SMILES code. Ignoring.")

# Definitions
def extract_ligand(pdb_file, ligand_name, output_file):
    '''Extracts the ligand from a MD PDB file and writes it to a new file.'''
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
    with open(output_file, 'w') as new_file:
        for line in lines:
            if ligand_name in line:
                new_file.write(line)

def remove_hydrogens_and_simplify(mol):
    '''Removes hydrogens, eliminates aromaticity and simplifies bonds in an RDKit molecule.'''
    mol = Chem.RemoveHs(mol)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.rdchem.BondType.SINGLE)
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    return mol

def mcs_match(smiles, ligand_pdb):
    '''Matches the Maximum Common Substructure (MCS) between a SMILES and a ligand PDB file.'''
    mol_from_pdb = Chem.MolFromPDBFile(ligand_pdb)
    mol_from_pdb = remove_hydrogens_and_simplify(mol_from_pdb)
    mol_from_smiles = Chem.MolFromSmiles(smiles)
    mol_from_smiles = remove_hydrogens_and_simplify(mol_from_smiles)
    mcs_result = rdFMCS.FindMCS([mol_from_pdb, mol_from_smiles])
    mcs_smarts = mcs_result.smartsString
    mol_from_mcs_smarts = Chem.MolFromSmarts(mcs_smarts)
    return mol_from_pdb.GetSubstructMatches(mol_from_mcs_smarts)  # List of lists with atom numbers

def return_atom_numbers(atom_matches, ligand_pdb):
    '''Returns the atom numbers in a ligand PDB file from atom matches.'''
    with open(ligand_pdb, 'r') as file:
        lines = file.readlines()
    atom_numbers = []
    for match in atom_matches[0]:  # Get the first list of atom numbers
        atom_numbers.append((int(lines[match].split()[1])) - 1)  # Get the atom number from the PDB file | -1 to match the 0-based indexing of MDTraj
    return atom_numbers

# Run the script
if smiles:
    action = 'inverse ' if inverse else ''
    print(f'This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files using the {action}maximum common substructure (MCS) with SMILES provided, with {reference_frame} being the reference frame.')
else:
    print(f'This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files, with {reference_frame} being the reference frame.')

# Initialize data structures for results
results_rmsd = pd.DataFrame(columns=['Ligand', 'Mean_RMSD', 'Std_RMSD'])
rmsd_data = defaultdict(list)

for subdir in os.listdir('.'):
    if not os.path.isdir(subdir):
        continue

    print(f'Analyzing {subdir}')

    pdb_file_path = os.path.join(subdir, pdb_file)
    xtc_file = os.path.join(subdir, traj_file)

    if not os.path.exists(pdb_file_path) or not os.path.exists(xtc_file):
        print(f'Trajectory or PDB file not found in {subdir}.')
        continue

    # Load trajectory excluding membrane and solvent to improve performance
    selection_query = 'not resname POPC and not resname SOL'
    selected_indices_traj = md.load(pdb_file_path).top.select(selection_query)
    traj = md.load(xtc_file, top=pdb_file_path, atom_indices=selected_indices_traj)

    # Select ligand atoms
    if smiles:
        ligand_pdb = f'{ligand_name}.pdb'
        extract_ligand(pdb_file_path, ligand_name, ligand_pdb)
        matching_atom_numbers = return_atom_numbers(mcs_match(smiles, ligand_pdb), ligand_pdb)
    else:
        matching_atom_numbers = traj.top.select(f'resname {ligand_name} and not element H')  # Exclude hydrogens

    if inverse and smiles:
        all_ligand_atoms = set(traj.top.select(f'resname {ligand_name} and not element H'))
        non_matching_atom_numbers = list(all_ligand_atoms - set(matching_atom_numbers))
        if len(non_matching_atom_numbers) == 0:
            print(f'No non-matching atoms found for ligand {ligand_name} in {subdir} using the provided SMILES with inverse flag.')
            continue
        lig_atom_numbers = non_matching_atom_numbers
    else:
        lig_atom_numbers = matching_atom_numbers

    if len(lig_atom_numbers) == 0:
        print(f'No atoms found for ligand {ligand_name} in {subdir}.')
        continue
    
    # Calculate RMSD
    traj = traj.atom_slice(lig_atom_numbers)
    rmsd_values = md.rmsd(traj, traj, reference_frame) * 10  # Convert to Angstrom
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    results_rmsd = results_rmsd.append({'Ligand': subdir, 'Mean_RMSD': rmsd_mean, 'Std_RMSD': rmsd_std}, ignore_index=True)
    print(f'Mean RMSD: {rmsd_mean:.2f}, Std RMSD: {rmsd_std:.2f}')

# Save results
results_rmsd = results_rmsd.sort_values(by='Ligand')
results_rmsd.to_csv(output_filename_rmsd, sep='\t', index=False)
print(f'RMSD results saved to {output_filename_rmsd}.')