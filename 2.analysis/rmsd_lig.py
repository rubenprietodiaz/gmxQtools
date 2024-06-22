import os
import argparse
import mdtraj as md # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from collections import defaultdict

# Argument parser setup and input parameters
parser = argparse.ArgumentParser(description="Calculate RMSD of a ligand in trajectory files.")
parser.add_argument('-p', '--pdb_file', type=str, default='finalOutput/start.pdb', help='Path to the PDB file. Default is finalOutput/start.pdb')
parser.add_argument('-l', '--ligand_name', type=str, default='L01', help='Name of the ligand. Default is L01')
parser.add_argument('-t', '--traj_file', type=str, default='traj_prod.xtc', help='Path to the trajectory file. Default is traj_prod.xtc')
parser.add_argument('-o', '--output_filename_rmsd', type=str, default='rmsd_lig_stat.txt', help='Output filename for RMSD results. Default is rmsd_stat.txt')
parser.add_argument('-f', '--reference_frame', type=int, default=0, help='Frame to use as reference for RMSD calculation. Default is 0')

args = parser.parse_args()

pdb_file = args.pdb_file
ligand_name = args.ligand_name
traj_file = args.traj_file
output_filename_rmsd = args.output_filename_rmsd
reference_frame = args.reference_frame

print(f"This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files using frame {reference_frame} as reference.")

results_rmsd = pd.DataFrame(columns=['Ligand', 'Mean_RMSD', 'Std_RMSD'])
rmsd_data = defaultdict(list) # Dictionary to store RMSD data for each group

for subdir in os.listdir('.'):
    if not os.path.isdir(subdir):
        continue

    print(f'Analyzing {subdir}...')

    pdb_file_path = os.path.join(subdir, pdb_file)
    xtc_file = os.path.join(subdir, traj_file)
    
    if not os.path.exists(pdb_file_path) or not os.path.exists(xtc_file): # If the files are not found, skip the directory
        print(f'Trajectory or PDB file not found in {subdir}.')
        continue

    # Load trajectory excluding membrane and solvent to improve performance
    selection_query = 'not resname POPC and not resname SOL'
    selected_indices_traj = md.load(pdb_file_path).top.select(selection_query)
    traj = md.load(xtc_file, top=pdb_file_path, atom_indices=selected_indices_traj)

    # Select ligand atoms
    ligand_atom_indices = traj.top.select(f'resname {ligand_name}')
    if len(ligand_atom_indices) == 0:
        print(f'No atoms found for ligand {ligand_name} in {subdir}.')
        continue

    # Calculations
    traj_ligand = traj.atom_slice(ligand_atom_indices)
    rmsd_values = md.rmsd(traj_ligand, traj_ligand, reference_frame) * 10 # To convert to Angstrom (from nm)
    rmsd_data[subdir.split('_')[0]].append(rmsd_values)
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    results_rmsd = results_rmsd.append({'Ligand': subdir, 'Mean_RMSD': rmsd_mean, 'Std_RMSD': rmsd_std}, ignore_index=True)
    print(f'Mean RMSD: {rmsd_mean:.2f}, Std RMSD: {rmsd_std:.2f}') # Print results with 2 decimal places for cleaning terminal output

# Save RMSD data for each group to a single .xvg file
for assay, rmsd_lists in rmsd_data.items():
    xvg_filename_rmsd = f"{assay}_rmsd.xvg"
    with open(xvg_filename_rmsd, 'w') as f:
        f.write('@ s0 legend "RMSD of {} for {}"\n'.format(ligand_name, assay))
        f.write('@ title "RMSD over time for {}"\n'.format(assay))
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
