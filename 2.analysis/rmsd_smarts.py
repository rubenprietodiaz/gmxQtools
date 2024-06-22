import os
import argparse
import mdtraj as md # type: ignore
import pandas as pd # type: ignore
import numpy as np  # type: ignore
from rdkit import Chem # type: ignore
from rdkit.Chem import AllChem # type: ignore
from collections import defaultdict

# Argument parser setup and input parameters
parser = argparse.ArgumentParser(description="Calculate RMSD of a ligand in trajectory files.")
parser.add_argument('-p', '--pdb_file', type=str, default='pymol/start.pdb', help='Path to the PDB file. Default is pymol/start.pdb')
parser.add_argument('-s', '--smarts_pattern', type=str, help='SMARTS pattern. If not provided, analyze the whole ligand')
parser.add_argument('-l', '--ligand_name', type=str, default='L01', help='Name of the ligand. Default is L01')
parser.add_argument('-t', '--traj_file', type=str, default='pymol/traj_prod_pymol.xtc', help='Path to the trajectory file. Default is pymol/traj_prod_pymol.xtc')
parser.add_argument('-o', '--output_filename_rmsd', type=str, default='rmsd_stat.txt', help='Output filename for RMSD results. Default is rmsd_stat.txt')
parser.add_argument('-f', '--reference_frame', type=int, default=0, help='Frame to use as reference for RMSD calculation. Default is 0')

args = parser.parse_args()

pdb_file = args.pdb_file
smarts_pattern = args.smarts_pattern
ligand_name = args.ligand_name
traj_file = args.traj_file
output_filename_rmsd = args.output_filename_rmsd
reference_frame = args.reference_frame

# Definitions
def extract_ligand(pdb_file, ligand_name, output_file):
    '''Extracts the ligand from a MD PDB file and writes it to a new file.'''
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
    with open(output_file, 'w') as new_file:
        for line in lines:
            if ligand_name in line:
                new_file.write(line)

def smarts_match(smarts, ligand_pdb):
    '''Matches a SMARTS pattern to a ligand PDB file.'''
    mol_from_pdb = Chem.MolFromPDBFile(ligand_pdb)
    mol_from_smarts = Chem.MolFromSmarts(smarts)
    return mol_from_pdb.GetSubstructMatches(mol_from_smarts) # List of lists with atom numbers

def return_atom_numbers(smarts, ligand_pdb):
    '''Returns the atom numbers in a ligand PDB file that match the SMARTS pattern.'''
    rdkit_numbers = smarts_match(smarts, ligand_pdb)[0] # Get the first list of atom numbers
    with open(ligand_pdb, 'r') as file:
        lines = file.readlines()
    atom_numbers = []
    for number in rdkit_numbers:
        atom_numbers.append((int(lines[number].split()[1])) - 1) # Get the atom number from the PDB file | -1 to match the 0-based indexing of MDTraj
    return atom_numbers

# Run the script
if smarts_pattern:
    print(f'This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files using the SMARTS pattern {smarts_pattern}.')
else:
    print(f'This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files.')

results_rmsd = pd.DataFrame(columns=['Ligand', 'Mean_RMSD', 'Std_RMSD'])
rmsd_data = defaultdict(list) # Dictionary to store RMSD data for each group

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
    selection_query = 'not resname POPC and not resname SOL' # Better performance with this query
    selected_indices_traj = md.load(pdb_file_path).top.select(selection_query)
    traj = md.load(xtc_file, top=pdb_file_path, atom_indices=selected_indices_traj)

    # Select ligand atoms
    if smarts_pattern:
        ligand_pdb = f'{ligand_name}.pdb'
        extract_ligand(pdb_file_path, ligand_name, ligand_pdb) # Extract ligand from the trajectory
        lig_atom_numbers = return_atom_numbers(smarts_pattern, ligand_pdb) # List of atom numbers in the ligand PDB file that match the SMARTS pattern
    else:
        lig_atom_numbers = traj.top.select(f'resname {ligand_name}')

    if len(lig_atom_numbers) == 0:
        print(f'No atoms found for ligand {ligand_name} in {subdir}.')
        continue

    # Calculations
    traj = traj.atom_slice(lig_atom_numbers) # Slice the trajectory to include only the ligand atoms
    rmsd_values = md.rmsd(traj, traj, reference_frame) * 10 # Convert to Angstrom
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    results_rmsd = results_rmsd.append({'Ligand': subdir, 'Mean_RMSD': rmsd_mean, 'Std_RMSD': rmsd_std}, ignore_index=True)
    print(f'Mean RMSD: {rmsd_mean:.2f}, Std RMSD: {rmsd_std:.2f}')
    
    # Store RMSD data for each group (assuming subdir names are like 'group_replica')
    group_name = subdir.split('_')[0]
    rmsd_data[group_name].append(rmsd_values)

    xvg_filename_rmsd = os.path.join(subdir, f"{subdir}_rmsd.xvg")
    with open(xvg_filename_rmsd, 'w') as f:
        f.write('@ s0 legend "RMSD of {} in {}"\n'.format(ligand_name, subdir))
        f.write('@ title "RMSD over time"\n')
        f.write('@ xaxis label "Frame"\n')
        f.write('@ yaxis label "RMSD (Å)"\n')
        for frame, rmsd in enumerate(rmsd_values):
            f.write(f"{frame} {rmsd}\n")
    print(f'RMSD values saved to {xvg_filename_rmsd}.')

# Save RMSD data for each group to a single .xvg file
for group, rmsd_lists in rmsd_data.items():
    xvg_filename_rmsd = f"{group}_rmsd.xvg"
    with open(xvg_filename_rmsd, 'w') as f:
        f.write('@ s0 legend "RMSD of {} for {}"\n'.format(ligand_name, group))
        f.write('@ title "RMSD over time for {}"\n'.format(group))
        f.write('@ xaxis label "Frame"\n')
        f.write('@ yaxis label "RMSD (Å)"\n')
        max_frames = max(len(rmsd) for rmsd in rmsd_lists) # Get the maximum number of frames
        for frame in range(max_frames): # Iterate over frames
            line = [f"{frame}"]
            for rmsd in rmsd_lists:
                if frame < len(rmsd):
                    line.append(f"{rmsd[frame]}")
                else:
                    line.append("")  # Fill missing values with empty strings
            f.write(" ".join(line) + "\n")

# To do: Add calculations for statistics of RMSD combined values (between replicas)

# Save results to a file
results_rmsd = results_rmsd.sort_values(by='Ligand')
results_rmsd.to_csv(output_filename_rmsd, sep='\t', index=False)
print(f'RMSD results saved to {output_filename_rmsd}.')
