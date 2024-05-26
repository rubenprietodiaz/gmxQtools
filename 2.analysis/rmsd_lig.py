import os
import mdtraj as md
import pandas as pd
import numpy as np

# Input
pdb_file = 'start.pdb'
traj_file = 'traj_prod.xtc'
output_filename_rmsd = 'rmsd_lig.txt'
output_filename_rmsf = 'rmsf_lig.txt'

# Definitions
results_rmsd = pd.DataFrame(columns=['Ligand', 'Mean_RMSD', 'Std_RMSD'])

for subdir in os.listdir('.'):
    if not os.path.isdir(subdir):
        continue

    print(f'Entering {subdir}...')

    # File paths setup
    pdb_file_path = os.path.join(subdir, pdb_file)
    xtc_file = os.path.join(subdir, traj_file)
    print(f'PDB file: {pdb_file_path}\nTrajectory file: {xtc_file}')

    # Files existence check
    if not os.path.exists(pdb_file_path) or not os.path.exists(xtc_file):
        print(f'Trajectory or PDB file not found in {subdir}.')
        continue

    # Load trajectory excluding unnecessary residues
    traj = md.load(xtc_file, top=pdb_file_path)
    print('Trajectory loaded successfully.')

    # Calculate RMSD
    rmsd_values = md.rmsd(traj, traj, 0) * 10  # Convert to Angstrom
    print(f'RMSD values: {rmsd_values}')
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    results_rmsd = results_rmsd.append({'Ligand': subdir, 'Mean_RMSD': rmsd_mean, 'Std_RMSD': rmsd_std}, ignore_index=True)
    print(f'Mean RMSD: {rmsd_mean}, Std RMSD: {rmsd_std}')

    # Save RMSD values to file
    xvg_filename_rmsd = os.path.join(subdir, f"{subdir}_rmsd.xvg")
    with open(xvg_filename_rmsd, 'w') as f:
        f.write('@ s0 legend "RMSD of {} in {}"\n'.format(subdir, subdir))
        f.write('@ title "RMSD over time"\n')
        f.write('@ xaxis label "Frame"\n')
        f.write('@ yaxis label "RMSD (Å)"\n')
        for frame, rmsd in enumerate(rmsd_values):
            f.write(f"{frame} {rmsd}\n")
    print(f'RMSD values saved to {xvg_filename_rmsd}.')


# Save results to files
results_rmsd = results_rmsd.sort_values(by='Ligand')
results_rmsd.to_csv(output_filename_rmsd, sep='\t', index=False)
print(f'RMSD results saved to {output_filename_rmsd}.')
