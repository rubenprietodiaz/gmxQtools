## UNDER CONSTRUCTION ##

import os
import shutil
import argparse
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description="Setup your md simulation after pymemdyn equilibration.")
    parser.add_argument("-d", "--dir", help="Directory of the complex to prepare FEP files. Omit if you don't want to prepare FEP files.")
    return parser.parse_args()

def prepare_fep_files(complex_directory, destination_folder_fep):
    """Prepare the FEP files for the complex in the specified directory."""
    # Check if the complex directory exists
    complex_directory = os.path.abspath(complex_directory)
    if not os.path.isdir(complex_directory):
        print(f"Error: The specified complex directory {complex_directory} does not exist.")
        return

    # Check if the destination folder for FEP files exists
    eqProd_directory = os.path.join(complex_directory, "eqProd")
    confout_gro = os.path.join(eqProd_directory, "confout200.gro")

    # Check if confout200.gro exists
    if not os.path.isfile(confout_gro):
        print(f"Error: confout200.gro not found in {eqProd_directory}.")
        return

    # Copy confout200.gro (before full relax) to ../confout_fep.gro
    confout_fep_gro = os.path.join(complex_directory, "confout_fep.gro")
    shutil.copy(confout_gro, confout_fep_gro)

    # Run gmx trjconv 
    command = ["gmx", "trjconv", "-pbc", "mol", "-center", "-ur", "compact", "-f", confout_fep_gro, "-o", "confout_fep.pdb"]
    with open(os.path.join(destination_folder_fep, "trjconv.log"), "w") as logfile:
        subprocess.run(command, stdout=logfile, stderr=subprocess.STDOUT, text=True, input="1\n0\n", cwd=complex_directory)

    # Move confout_fep.pdb to destination_folder_fep
    shutil.move(os.path.join(complex_directory, "confout_fep.pdb"), os.path.join(destination_folder_fep, "confout_fep.pdb"))
    print(f"File confout_fep.pdb created successfully in {destination_folder_fep}.")

def process_water_molecules(input_pdb, output_pdb):
    """Process water molecules in the input PDB file. Generate a new PDB file with only water molecules around the polymer."""
    pymol_script = f"""
load {input_pdb}
select water, solvent within 10 of polymer
select water_2, solvent within 10 of resn popc
select water_3, water and water_2
save {output_pdb}, water_3
"""
    with open("get_waters.pml", "w") as file:
        file.write(pymol_script)
    subprocess.run(["pymol", "-c", "get_waters.pml"], check=True)
    os.remove("get_waters.pml")

def rename_and_remove_elements(pdb_file):
    """Rename and remove elements in the PDB file."""
    subprocess.run(["grep", "'SOL'", pdb_file, ">", "temp.pdb"])
    shutil.move("temp.pdb", pdb_file.replace("confout_water", "confout"))
    subprocess.run(["sed", "-i", "'s/SOL/HOH/g'", pdb_file])
    subprocess.run(["sed", "-i", "'s/OW/O /g'", pdb_file])
    subprocess.run(["sed", "-i", "'s/HW1/H1 /g'", pdb_file])
    subprocess.run(["sed", "-i", "'s/HW2/H2 /g'", pdb_file])
    subprocess.run(["sed", "-i", "'s/O2 /OXT/g'", pdb_file]) # Check this modifications
    

def clean_pdb(original_pdb, ligand):
    temp_file = original_pdb.replace('.pdb', '_temp.pdb')
    ligand_file = os.path.join(os.path.dirname(original_pdb), f"{ligand}.pdb")

    with open(original_pdb, 'r') as file:
        lines = file.readlines()

    with open(temp_file, 'w') as out_file, open(ligand_file, 'w') as ligand_out:
        for line in lines:
            if ligand in line:
                ligand_out.write(line)
            if ligand not in line and 'SOL' not in line and 'CL-' not in line and 'O1 ' not in line and 'END' not in line and 'TER' not in line:
                out_file.write(line)
    
    shutil.move(temp_file, original_pdb)
    print(f"Removed unwanted lines from {original_pdb} and saved ligand to {ligand_file}.")
    
def merge_water_and_pdb(pdb_file, water_pdb):
    """Merge the water molecules and the PDB file."""
    with open(water_pdb, 'r') as water_file, open(pdb_file, 'a') as file:
        file.writelines(water_file.readlines())

def remove_last_two_characters(original_file, modified_file):
    """Remove the last two characters from each line in the original file and save the modified file, for better visualization."""
    with open(original_file, 'r') as file_in, open(modified_file, 'w') as file_out:
        for line in file_in:
            if line.strip():
                file_out.write(line[:-2] + '\n')

# RUN
args = parse_arguments()

destination_folder_fep = "fep_preparation_files"
os.makedirs(destination_folder_fep, exist_ok=True)
prepare_fep_files(args.dir, destination_folder_fep)

# Define paths to the PDB files
confout_pdb = os.path.join(destination_folder_fep, "confout_fep.pdb")
water_pdb = os.path.join(destination_folder_fep, "confout_water.pdb")

# Process water molecules
process_water_molecules(confout_pdb, water_pdb)
clean_pdb(confout_pdb, "L01")
merge_water_and_pdb(confout_pdb, water_pdb)

os.remove(water_pdb) # Remove the temporary file of waters