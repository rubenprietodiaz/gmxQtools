import os
import sys
import shutil
import fileinput
import argparse
import subprocess
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Setup your md simulation after pymemdyn equilibration.")
    parser.add_argument("-d", "--dir", help="Directory of the complex to prepare FEP files. Omit if you don't want to prepare FEP files.")
    parser.add_argument("-nc", "--noclean", help="Do not remove logs and temporary files.", action="store_true")
    return parser.parse_args()

def prepare_fep_files(complex_directory, destination_folder_fep):
    """Generate and organize the necessary files for FEP calculations."""
    complex_directory = os.path.abspath(complex_directory)

    if not os.path.isdir(complex_directory):
        print(f"ERROR: The specified complex directory {complex_directory} does not exist. Possible problem with the name of input files.")
        return
    eqProd_directory = os.path.join(complex_directory, "eqProd")
    confout_gro = os.path.join(eqProd_directory, "confout200.gro")

    if not os.path.isfile(confout_gro):
        print(f"ERROR: confout200.gro not found in {eqProd_directory}. Possible error running the equilibration. Check the log files.")
        return

    # Copy confout200.gro (before full relax) to ../confout_fep.gro
    confout_fep_gro = os.path.join(complex_directory, "confout_fep.gro")
    shutil.copy(confout_gro, confout_fep_gro)

    # Run gmx trjconv 
    command = ["gmx", "trjconv", "-pbc", "mol", "-center", "-ur", "compact", "-f", confout_fep_gro, "-o", "tmp.pdb"] # Generate temporary pdb file from confout_fep.gro
    with open(os.path.join(destination_folder_fep, "trjconv.log"), "w") as logfile: # Save trjconv output to logfile
        subprocess.run(command, stdout=logfile, stderr=subprocess.STDOUT, text=True, input="1\n0\n", cwd=complex_directory)

    # Move temporary file to destination_folder_fep
    shutil.move(os.path.join(complex_directory, "tmp.pdb"), os.path.join(destination_folder_fep, "tmp.pdb"))
    
def remove_last_two_characters(original_file, modified_file):
    """Remove the last two characters of pdb files for avoid errors."""
    with open(original_file, 'r') as file_in, open(modified_file, 'w') as file_out:
        for line in file_in:
            if line.strip():
                file_out.write(line[:-2] + '\n')

def process_complex(input_pdb): # 30/04/24 Removed select water, resn SOL and (polymer around 10) and (lig around 25) and divided into two selections
    """Process water molecules in the input PDB file. Generate a new PDB file with only water molecules around the polymer."""
    pymol_script = f"""
load {input_pdb}
select lig, resn L01 
select water, resn SOL and (polymer around 10)
select water2, resn SOL and (lig around 25)
select membrane, resn POPC and byres (lig around 32)
create complex, polymer or membrane
create water, water or water2
save {destination_folder_fep}/complex.pdb, complex
save {destination_folder_fep}/ligand.pdb, lig  
save {destination_folder_fep}/water.pdb, water
"""
    with open("process.pml", "w") as file:
        file.write(pymol_script)
    
    log_file_path = os.path.join(destination_folder_fep, "pymol.log")
    
    with open(log_file_path, "w") as log_file:
        subprocess.run(["pymol", "-c", "process.pml"], stdout=log_file, stderr=subprocess.STDOUT, check=True)
    os.remove("process.pml")
    os.remove(input_pdb)


# The following functions were taken and adapted from the script process_solvent.py, merge.py and fix_pdb_Q.sh 
# in the MD_analysis repository (https://github.com/bananatana/MD_analysis/tree/main/GROMACS_to_Q)

def read_pdb(filename):
    """Read the atoms from a PDB file and return a list of tuples with the atom type and coordinates."""
    atoms = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_type = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                if atom_type == 'OW': atom_type = 'O'
                if atom_type == 'HW': atom_type = 'H'
                if atom_type == 'HW1': atom_type = 'H'
                if atom_type == 'HW2': atom_type = 'H'
                atoms.append((atom_type, np.array([x, y, z])))
    return atoms

def compute_distances(atoms):
    """Compute the distances between oxygen and hydrogen atoms to identify water molecules."""
    oxygens = [atom for atom in atoms if atom[0] == 'O']
    hydrogens = [atom for atom in atoms if atom[0] == 'H']
    
    water_molecules = []
    for o_atom in oxygens:
        distances = [(h_atom, np.linalg.norm(o_atom[1] - h_atom[1])) for h_atom in hydrogens]
        distances.sort(key=lambda x: x[1])
        
        # Assuming the two closest hydrogens form a water molecule with the oxygen
        # and the distance is within a reasonable range for O-H bonds
        if distances[0][1] < 1.5 and distances[1][1] < 1.5:
            water_molecules.append((o_atom, distances[0][0], distances[1][0]))
    
    return water_molecules

def write_pdb(water_molecules, output_file):
    """Write the water molecules to a PDB file."""
    with open(output_file, 'w') as file:
        atom_id = 1
        res_id = 1
        for water in water_molecules:
            o_atom, h1_atom, h2_atom = water

            # Write oxygen atom
            file.write(f"ATOM  {atom_id:>5}  OW  HOH A{res_id:>4}    {o_atom[1][0]:>8.3f}{o_atom[1][1]:>8.3f}{o_atom[1][2]:>8.3f}  1.00  0.00           O\n")
            atom_id += 1

            # Write first hydrogen atom
            file.write(f"ATOM  {atom_id:>5} 1H   HOH A{res_id:>4}    {h1_atom[1][0]:>8.3f}{h1_atom[1][1]:>8.3f}{h1_atom[1][2]:>8.3f}  1.00  0.00           H\n")
            atom_id += 1

            # Write second hydrogen atom
            file.write(f"ATOM  {atom_id:>5} 2H   HOH A{res_id:>4}    {h2_atom[1][0]:>8.3f}{h2_atom[1][1]:>8.3f}{h2_atom[1][2]:>8.3f}  1.00  0.00           H\n")
            atom_id += 1

            res_id += 1

def merge_pdb_files(file1, file2, output):
    """Merge the content of two PDB files into a single file."""
    with open(file1, 'r') as f1, open(file2, 'r') as f2, open(output, 'w') as out:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

        # Write content of the first file
        for line in lines1[:-1]:  # Exclude the last line (END) of the first file
            out.write(line)

        # If there are no lines starting with "END" in the first file, add one
        if not any(line.startswith('END') for line in lines1):
            out.write('END\n')

        # Write content of the second file
        for line in lines2:
            out.write(line)

        # If there are no lines starting with "END" in the second file, add one
        if not any(line.startswith('END') for line in lines2):
            out.write('END\n')

def correct_pdb(file_path):
    substitutions = [
        ('POPC', 'POP '),
        ('1H   HOH', ' H1  HOH'),
        ('2H   HOH', ' H2  HOH'),
        ('OW  HOH', 'O   HOH'),
        ('CD  ILE', 'CD1 ILE'),
        ('O1 ', 'O  '), # For C-terminus
        ('O2 ', 'OXT'), # For C-terminus
        ('HISH', 'HIP '),
        ('HISD', 'HID '),
        ('HISE', 'HIE ')
    ]

    # Open the file for in-place editing
    with fileinput.input(files=[file_path], inplace=True, backup='.bak') as file:
        for line in file:
            for old, new in substitutions:
                line = line.replace(old, new)
            sys.stdout.write(line)  # Write directly to the file

def clean_up(directory):
    """Remove unnecessary files."""
    files_to_remove = ["confout_fep.gro", "trjconv.log", "tmp.pdb", "pymol.log", "water_fixed.pdb", "water.pdb", "system.pdb.bak", "complex.pdb"]
    for file in files_to_remove:
        file_path = os.path.join(directory, file)
        if os.path.isfile(file_path):
            os.remove(file_path)

# Running
args = parse_arguments()
print(f"Preparing FEP files for {args.dir}.")
destination_folder_fep = "fep_preparation_files"
os.makedirs(destination_folder_fep, exist_ok=True)
print(f"[1/11]  Directory created: {destination_folder_fep}.")
prepare_fep_files(args.dir, destination_folder_fep)
print(f"[2/11]  Necessary files from {args.dir} copied to {destination_folder_fep}.")
print(f"[3/11]  PDB files generated.")
remove_last_two_characters(os.path.join(destination_folder_fep, "tmp.pdb"), os.path.join(destination_folder_fep, "tmp_clean.pdb"))
print(f"[4/11]  Last two characters removed from PDB files.")
shutil.move(os.path.join(destination_folder_fep, "tmp_clean.pdb"), os.path.join(destination_folder_fep, "tmp.pdb"))
process_complex(os.path.join(destination_folder_fep, "tmp.pdb"))
print(f"[5/11]  Complex processed: complex (protein & membrane), waters and ligand saved.")
atoms = read_pdb(os.path.join(destination_folder_fep, "water.pdb"))
print(f"[6/11]  Water molecules read. Starting re-numbering. It may take long time.")
water_molecules = compute_distances(atoms)
print(f"[7/11]  Water molecules numbering re-calculated.")
write_pdb(water_molecules, os.path.join(destination_folder_fep, 'water_fixed.pdb'))
print(f"[8/11]  Water molecules written to water_fixed.pdb.")
merge_pdb_files(os.path.join(destination_folder_fep, 'complex.pdb'), os.path.join(destination_folder_fep, 'water_fixed.pdb'), os.path.join(destination_folder_fep, 'system.pdb'))
print(f"[9/11]  Protein, membrane and fixed waters merged into system.pdb.")
correct_pdb(os.path.join(destination_folder_fep, 'system.pdb'))
print(f"[10/11] PDB file corrected.")

if args.noclean:
    print(f"[11/11] Temporary files, logs and scripts not removed. FEP files prepared in {destination_folder_fep}.")
else:
    clean_up(destination_folder_fep)
    print(f"[11/11] Temporary files, logs and scripts removed. FEP files prepared in {destination_folder_fep}.\nRelevant files: system.pdb and ligand.pdb.")