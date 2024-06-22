import os
import mdtraj as md # type: ignore
import pandas as pd # type: ignore
import numpy as np  # type: ignore
from rdkit import Chem # type: ignore
from rdkit.Chem import AllChem # type: ignore

# Input // change to argparse (add change smarts-smiles? or maximmum common substructure?)
pdb_file = 'pymol/start.pdb'
ref_pdb_file = 'ref.pdb'
smarts_pattern = '[#8]-[#6]-[#7]-[#6]1-[#6](-[#6]-[#6](-[#6]-[#6]-1-[#9])-[#6](-[#8])-[#7]-[#6]1-[#7]-[#6]-[#6]-[#16]-1)-[#9]' # for Lundbeck project
ligand_name = 'L01'
ligand_pdb = f'{ligand_name}.pdb'
ref_pdb = 'ref_lig.pdb'
traj_file = 'pymol/traj_prod_pymol.xtc'
xvg_filename_rmsd = 'rmsd.xvg'
xvg_filename_rmsf = 'rmsf.xvg'
output_filename_rmsd = 'rmsd_results_old.txt'
output_filename_rmsf = 'rmsf_results_old.txt'

# To do: make possible to input the SMARTS pattern with H (remove them) and multiple bonds (e.g. aromatic bonds)

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

def calculations(value_type, traj, ref_traj): # Unused
    '''Calculates different values for a ligand in a trajectory. Value types: rmsd, rmsf.'''
    if value_type == 'rmsd':
        values = md.rmsd(traj, ref_traj, 0) * 10
        print(f'RMSD values: {values}')

    elif value_type == 'rmsf':
        values = md.rmsf(traj, ref_traj, 0) * 10
        print(f'RMSF values: {values}')
    return values

def statistics(value_type, traj, ref_traj): # Unused
    '''Calculates the mean and standard deviation of RMSD or RMSF values.'''
    values = calculations(value_type, traj, ref_traj)
    mean = np.mean(values)
    std = np.std(values)
    print(f'Mean: {mean}, Std: {std}')
    if value_type == 'rmsd':
        results_rmsd = results_rmsd.append({'Ligand': subdir, f'Mean_RMSD': mean, f'Std_RMSD': std}, ignore_index=True)
    elif value_type == 'rmsf':
        results_rmsf = results_rmsf.append({'Ligand': subdir, f'Mean_RMSF': mean, f'Std_RMSF': std}, ignore_index=True)
    return results_rmsd, results_rmsf

results_rmsd = pd.DataFrame(columns=['Ligand', 'Mean_RMSD', 'Std_RMSD'])
results_rmsf = pd.DataFrame(columns=['Ligand', 'Mean_RMSF', 'Std_RMSF'])

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
    selection_query = 'not resname POPC and not resname SOL' # Better performance with this query
    selected_indices_traj = md.load(pdb_file_path).top.select(selection_query)
    # selected_indices_ref = md.load(ref_pdb_file).top.select(selection_query)
    traj = md.load(xtc_file, top=pdb_file_path, atom_indices=selected_indices_traj)
    # ref_traj = md.load(ref_pdb_file, atom_indices=selected_indices_ref)
    print('Trajectory and reference pose loaded successfully.')

    # Process the ligand
    extract_ligand(pdb_file_path, ligand_name, ligand_pdb) # Extract ligand from the trajectory
    # extract_ligand(ref_pdb_file, ligand_name, ref_pdb) # Extract ligand from the reference pose

    lig_atom_numbers = return_atom_numbers(smarts_pattern, ligand_pdb) # List of atom numbers in the ligand PDB file that match the SMARTS pattern
    # ref_atom_numbers = return_atom_numbers(smarts_pattern, ref_pdb) # List of atom numbers in the ligand PDB file that match the SMARTS pattern

    # Align the trajectory with the reference
    # alignment_indices = traj.top.select('backbone')
    # traj.superpose(ref_traj, atom_indices=alignment_indices)
    # print('Trajectory aligned with the reference PDB.')

    # Slice the trajectory to include only the ligand coinciding SMARTS pattern
    traj = traj.atom_slice(lig_atom_numbers)
    # ref_traj = ref_traj.atom_slice(ref_atom_numbers)

    # Calculate RMSD
    rmsd_values = md.rmsd(traj, traj, 0) * 10 # Convert to Angstrom || If ref_traj is used, the RMSD is calculated between the two trajectories
    print(f'RMSD values: {rmsd_values}')
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    results_rmsd = results_rmsd.append({'Ligand': subdir, 'Mean_RMSD': rmsd_mean, 'Std_RMSD': rmsd_std}, ignore_index=True)
    print(f'Mean RMSD: {rmsd_mean}, Std RMSD: {rmsd_std}')

    xvg_filename_rmsd = os.path.join(subdir, f"{subdir}_rmsd.xvg")
    with open(xvg_filename_rmsd, 'w') as f:
        f.write('@ s0 legend "RMSD of {} in {}"\n'.format(ligand_name, subdir))
        f.write('@ title "RMSD over time"\n')
        f.write('@ xaxis label "Frame"\n')
        f.write('@ yaxis label "RMSD (Å)"\n')
        for frame, rmsd in enumerate(rmsd_values):
            f.write(f"{frame} {rmsd}\n")
    print(f'RMSD values saved to {xvg_filename_rmsd}.')

    # Calculate RMSF
    rmsf_values = md.rmsf(traj, traj, 0) * 10 # Convert to Angstrom
    print(f'RMSF values: {rmsf_values}')
    rmsf_mean = np.mean(rmsf_values)
    rmsf_std = np.std(rmsf_values)
    results_rmsf = results_rmsf.append({'Ligand': subdir, 'Mean_RMSF': rmsf_mean, 'Std_RMSF': rmsf_std}, ignore_index=True)
    print(f'Mean RMSF: {rmsf_mean}, Std RMSF: {rmsf_std}')

    xvg_filename_rmsf = os.path.join(subdir, f"{subdir}_rmsf.xvg")
    with open(xvg_filename_rmsf, 'w') as f:
        f.write('@ s0 legend "RMSF of {} in {}"\n'.format(ligand_name, subdir))
        f.write('@ title "RMSF over time"\n')
        f.write('@ xaxis label "Residue"\n')
        f.write('@ yaxis label "RMSF (Å)"\n')
        for residue, rmsf in enumerate(rmsf_values):
            f.write(f"{residue} {rmsf}\n")
    print(f'RMSF values saved to {xvg_filename_rmsf}.')

    # Print selected atoms for rmsd calculation (sanity check). Decreases performance.
    # To do: make this optional instead of always printing
    # sliced_traj = traj.atom_slice(lig_atom_numbers)
    # for atom in sliced_traj.topology.atoms:
    #     print(f"Atom index: {atom.index}, Atom name: {atom.name}, Element: {atom.element.symbol}, Residue name: {atom.residue.name}, Residue index: {atom.residue.index}")

    # ref_sliced_traj = traj.atom_slice(ref_atom_numbers)
    # for atom in ref_sliced_traj.topology.atoms:
    #     print(f"Atom index: {atom.index}, Atom name: {atom.name}, Element: {atom.element.symbol}, Residue name: {atom.residue.name}, Residue index: {atom.residue.index}")

    print(f'Exiting {subdir}...')

# Save results to a file
results_rmsd = results_rmsd.sort_values(by='Ligand')
results_rmsd.to_csv(output_filename_rmsd, sep='\t', index=False)
print(f'RMSD results saved to {output_filename_rmsd}.')

results_rmsf = results_rmsf.sort_values(by='Ligand')
results_rmsf.to_csv(output_filename_rmsf, sep='\t', index=False)
print(f'RMSF results saved to {output_filename_rmsf}.')