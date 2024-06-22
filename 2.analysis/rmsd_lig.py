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
parser.add_argument('-p', '--pdb_file', type=str, default='pymol/start.pdb', help='Path to the PDB file. Default is pymol/start.pdb')  # Change to finalOutput/start.pdb
parser.add_argument('-s', '--smarts_pattern', type=str, help='SMARTS pattern. If not provided, analyze the whole ligand. When used, please provide a simplified SMARTS pattern (no Hs, no aromaticity, no bond orders)')
parser.add_argument('-S', '--smiles', type=str, help='SMILES code. If provided, calculate maximum common substructure (MCS) with the SMILES. If you want to analyze a specific part of the ligand, you can provide the SMILES code for that part.')
parser.add_argument('-i', '--inverse', action='store_true', help='Calculate RMSD for atoms not matching the SMARTS pattern or MCS with SMILES')
parser.add_argument('-l', '--ligand_name', type=str, default='L01', help='Name of the ligand. Default is L01')
parser.add_argument('-t', '--traj_file', type=str, default='pymol/traj_prod_pymol.xtc', help='Path to the trajectory file. Default is pymol/traj_prod_pymol.xtc')  # Change to finalOutput/traj_prod_pymol.xtc
parser.add_argument('-o', '--output_filename_rmsd', type=str, default='rmsd_stat.txt', help='Output filename for RMSD results. Default is rmsd_stat.txt')
parser.add_argument('-f', '--reference_frame', type=int, default=0, help='Frame to use as reference for RMSD calculation. Default is 0')

args = parser.parse_args()

pdb_file = args.pdb_file
smarts_pattern = args.smarts_pattern
smiles = args.smiles
inverse = args.inverse
ligand_name = args.ligand_name
traj_file = args.traj_file
output_filename_rmsd = args.output_filename_rmsd
reference_frame = args.reference_frame

# Check for conflicting arguments
if inverse and not smarts_pattern and not smiles:
    print("Inverse flag provided without SMARTS pattern or SMILES code. Ignoring.")

if smarts_pattern and smiles:
    raise ValueError("Cannot provide both SMARTS pattern and SMILES code. Please choose one.")

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

def smarts_match(smarts, ligand_pdb): # Depcrecated
    '''Matches a SMARTS pattern to a ligand PDB file.'''
    mol_from_pdb = Chem.MolFromPDBFile(ligand_pdb)
    mol_from_pdb = remove_hydrogens_and_simplify(mol_from_pdb)
    mol_from_smarts = Chem.MolFromSmarts(smarts)
    mol_from_smarts = remove_hydrogens_and_simplify(mol_from_smarts)
    return mol_from_pdb.GetSubstructMatches(mol_from_smarts)  # List of lists with atom numbers

def mcs_match(smiles, ligand_pdb): # Remove hydrogens and simplify for matching
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
if smarts_pattern:
    action = 'inverse ' if inverse else ''
    print(f'This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files using the {action}SMARTS pattern provided, with {reference_frame} being the reference frame.')
elif smiles:
    action = 'inverse ' if inverse else ''
    print(f'This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files using the {action}maximum common substructure (MCS) with SMILES provided, with {reference_frame} being the reference frame.')
else:
    print(f'This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files, with {reference_frame} being the reference frame.')

results_rmsd = pd.DataFrame(columns=['Ligand', 'Mean_RMSD', 'Std_RMSD'])
rmsd_data = defaultdict(lambda: defaultdict(list))  # Dictionary to store RMSD data for each group and ligand

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
    selection_query = 'not resname POPC and not resname SOL'  # Better performance with this query
    selected_indices_traj = md.load(pdb_file_path).top.select(selection_query)
    traj = md.load(xtc_file, top=pdb_file_path, atom_indices=selected_indices_traj)

    # Select ligand atoms
    if smarts_pattern:
        ligand_pdb = f'{ligand_name}.pdb'
        extract_ligand(pdb_file_path, ligand_name, ligand_pdb)  # Extract ligand from the trajectory
        matching_atom_numbers = return_atom_numbers(smarts_match(smarts_pattern, ligand_pdb), ligand_pdb)  # List of atom numbers in the ligand PDB file that match the SMARTS pattern
    elif smiles:
        ligand_pdb = f'{ligand_name}.pdb'
        extract_ligand(pdb_file_path, ligand_name, ligand_pdb)  # Extract ligand from the trajectory
        matching_atom_numbers = return_atom_numbers(mcs_match(smiles, ligand_pdb), ligand_pdb)  # List of atom numbers in the ligand PDB file that match the MCS with SMILES
    else:
        matching_atom_numbers = traj.top.select(f'resname {ligand_name} and not element H')  # Excluir hidrógenos si no se da smarts o smiles

    if inverse and (smarts_pattern or smiles):
        all_ligand_atoms = set(traj.top.select(f'resname {ligand_name} and not element H'))  # Excluir hidrógenos
        non_matching_atom_numbers = list(all_ligand_atoms - set(matching_atom_numbers))
        if len(non_matching_atom_numbers) == 0:
            print(f'No non-matching atoms found for ligand {ligand_name} in {subdir} using the provided SMARTS/SMILES with inverse flag.')
            continue
        lig_atom_numbers = non_matching_atom_numbers
    else:
        lig_atom_numbers = matching_atom_numbers

    if len(lig_atom_numbers) == 0:
        print(f'No atoms found for ligand {ligand_name} in {subdir}.')
        continue

    # Calculations
    traj = traj.atom_slice(lig_atom_numbers)  # Slice the trajectory to include only the ligand atoms
    rmsd_values = md.rmsd(traj, traj, reference_frame) * 10  # Convert to Angstrom
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    results_rmsd = results_rmsd.append({'Ligand': subdir, 'Mean_RMSD': rmsd_mean, 'Std_RMSD': rmsd_std}, ignore_index=True)
    print(f'Mean RMSD: {rmsd_mean:.2f}, Std RMSD: {rmsd_std:.2f}')
    
    # Store RMSD data for each group and ligand (assuming subdir names are like 'group_replica')
    group_name = subdir.split('_')[0]
    rmsd_data[group_name][ligand_name].append(rmsd_values)

    xvg_filename_rmsd = os.path.join(subdir, f"{subdir}_rmsd.xvg")
    with open(xvg_filename_rmsd, 'w') as f:
        f.write('@ s0 legend "RMSD of {} in {}"\n'.format(ligand_name, subdir))
        f.write('@ title "RMSD over time"\n')
        f.write('@ xaxis label "Frame"\n')
        f.write('@ yaxis label "RMSD (Å)"\n')
        for frame, rmsd in enumerate(rmsd_values):
            f.write(f"{frame} {rmsd}\n")
    print(f'RMSD values saved to {xvg_filename_rmsd}.')

# Calculate mean and standard deviation of RMSD for each frame across all replicas for each ligand
combined_rmsd_stats = []

for group, ligands in rmsd_data.items():
    max_frames = max(len(rmsd) for ligand in ligands.values() for rmsd in ligand)
    for frame in range(max_frames):
        row = {'Frame': frame}
        for ligand, rmsd_lists in ligands.items():
            frame_rmsd_values = [rmsd[frame] for rmsd in rmsd_lists if frame < len(rmsd)]
            if frame_rmsd_values:
                frame_mean = np.mean(frame_rmsd_values)
                frame_std = np.std(frame_rmsd_values)
                row[f'{ligand}_Mean_RMSD'] = frame_mean
                row[f'{ligand}_Std_RMSD'] = frame_std
                row[f'{ligand}_Replicas'] = len(frame_rmsd_values)
        combined_rmsd_stats.append(row)

combined_rmsd_df = pd.DataFrame(combined_rmsd_stats)
combined_rmsd_df.to_csv('combined_rmsd_stats_by_frame.txt', sep='\t', index=False)
print('Combined RMSD statistics by frame saved to combined_rmsd_stats_by_frame.txt.')

# Save results to a file
results_rmsd = results_rmsd.sort_values(by='Ligand')
results_rmsd.to_csv(output_filename_rmsd, sep='\t', index=False)
print(f'RMSD results saved to {output_filename_rmsd}.')
