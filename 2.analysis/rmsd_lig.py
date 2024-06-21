import os
import mdtraj as md # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from collections import defaultdict

# Input
pdb_file = 'pymol/start.pdb' # Now should be 'finalOutput/start.pdb'
ligand_name = 'L01'
traj_file = 'traj_prod.xtc'
output_filename_rmsd = 'rmsd_stat.txt'

print(f"This script will calculate the RMSD of the ligand {ligand_name} in the trajectory files.")

results_rmsd = pd.DataFrame(columns=['Ligand', 'Mean_RMSD', 'Std_RMSD'])
rmsd_data = defaultdict(list) # Dictionary to store RMSD data for each group

for subdir in os.listdir('.'):
    if not os.path.isdir(subdir):
        continue

    print(f'Entering {subdir}...')

    pdb_file_path = os.path.join(subdir, pdb_file)
    xtc_file = os.path.join(subdir, traj_file)
    print(f'PDB file: {pdb_file_path}\nTrajectory file: {xtc_file}')

    if not os.path.exists(pdb_file_path) or not os.path.exists(xtc_file): # If the files are not found, skip the directory
        print(f'Trajectory or PDB file not found in {subdir}.')
        continue

    # Load trajectory excluding membrane and solvent to improve performance
    selection_query = 'not resname POPC and not resname SOL'
    selected_indices_traj = md.load(pdb_file_path).top.select(selection_query)
    traj = md.load(xtc_file, top=pdb_file_path, atom_indices=selected_indices_traj)
    print('Trajectory loaded successfully.')

    # Select ligand atoms
    ligand_atom_indices = traj.top.select(f'resname {ligand_name}')
    if len(ligand_atom_indices) == 0:
        print(f'No atoms found for ligand {ligand_name} in {subdir}.')
        continue

    # Calculate ligand RMSD
    traj_ligand = traj.atom_slice(ligand_atom_indices)
    rmsd_values = md.rmsd(traj_ligand, traj_ligand, 0) * 10 # To convert to Angstrom (from nm)
    rmsd_data[subdir.split('_')[0]].append(rmsd_values)
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    results_rmsd = results_rmsd.append({'Ligand': subdir, 'Mean_RMSD': rmsd_mean, 'Std_RMSD': rmsd_std}, ignore_index=True)
    print(f'Mean RMSD: {rmsd_mean}, Std RMSD: {rmsd_std}')

    print(f'Exiting {subdir}...')

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
    print(f'Grouped RMSD values saved to {xvg_filename_rmsd}.')

# Save results to a file
results_rmsd = results_rmsd.sort_values(by='Ligand')
results_rmsd.to_csv(output_filename_rmsd, sep='\t', index=False)
print(f'RMSD results saved to {output_filename_rmsd}.')
