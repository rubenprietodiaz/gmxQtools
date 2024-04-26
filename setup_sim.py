import os
import shutil
import argparse
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to copy files, create submit.sh, and modify simulation duration.")
    parser.add_argument("-t", "--simulation-time", type=int, help="Simulation time in nanoseconds (default by PyMemDyn: 10 ns).")
    parser.add_argument("-rt", "--runtime", type=int, help="Runtime in hours (default: 24).", default=24)
    parser.add_argument("-C", "--cluster", choices=["CSB", "CESGA", "TETRA"], default="TETRA", help="Choose the cluster (default: TETRA).")
    parser.add_argument("-f", "--fep", help="Directory of the complex to prepare FEP files. Omit if you don't want to prepare FEP files.")
    return parser.parse_args()

def copy_files_in_directory(directory, destination_folder):
    """Iterate over the folders within the current directory and call the copy_files function."""
    # Iterate over the folders within the current directory
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            current_folder = os.path.join(root, dir)
            # Call the copy_files function to copy files in the current folder
            copy_files(current_folder, destination_folder)
    # Create submit script outside the loop
    create_submit_script(destination_folder)
    print(f"Folder {destination_folder} with scripts and necessary files created. Run 'cd {destination_folder}' and 'sh submit.sh' to start the simulations. If you want to do more than one replica, copy {destination_folder} before running any simulation.")


def copy_files(folder, destination_folder):
    """Copy the required files for MD simulations to the destination folder."""
    # Check the existence of each required file
    if all(os.path.isfile(os.path.join(folder, file)) for file in ["prod.mdp", "topol.top", "index.ndx","topol.tpr"]) and os.path.exists(os.path.join(folder, "finalOutput", "confout.gro")):
        # Extract the name of the current folder
        folder_name = os.path.basename(folder)

        # Create the folder for files within the destination folder
        destination_folder_path = os.path.join(destination_folder, folder_name)
        os.makedirs(destination_folder_path, exist_ok=True)

        # Copy the files to the corresponding folder
        for file in ["prod.mdp", "topol.top", "index.ndx", "topol.tpr"]:
            shutil.copy(os.path.join(folder, file), destination_folder_path)

        # Check the existence of the finalOutput folder
        final_output_folder = os.path.join(folder, "finalOutput")
        if os.path.isdir(final_output_folder):
            confout_gro = os.path.join(final_output_folder, "confout.gro")
            if os.path.isfile(confout_gro):
                shutil.copy(confout_gro, destination_folder_path)

        # Copy .itp files if they exist
        for itp_file in os.listdir(folder):
            if itp_file.endswith(".itp"):
                shutil.copy(os.path.join(folder, itp_file), destination_folder_path)

        # Create the run_md.sh script inside the destination folder
        create_run_md_script(destination_folder_path)

        # Modify prod.mdp if simulation time is provided
        if args.simulation_time:
            modify_simulation_time(destination_folder_path)
    else:
        return False

def create_submit_script(destination_folder):
    """Create the submit.sh script to submit the MD simulations."""
    submit_script_content = """#!/bin/bash

start_dir=$(pwd)

cd "$1"
for dir in */ ; do
    cd "$dir"
    sbatch run_md.sh
    cd "$start_dir/$1"
done
"""
    submit_script_path = os.path.join(destination_folder, "submit.sh")
    with open(submit_script_path, "w") as submit_script_file:
        submit_script_file.write(submit_script_content)

def create_run_md_script(destination_folder):
    """Create the run_md.sh script to run the MD simulations in all the subfolders."""
    if args.cluster == "CSB":
        run_md_script_content = f"""#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t {args.runtime}:00:00
#SBATCH --gpus-per-task=1
#SBATCH --job-name=pymemdyn
#              d-hh:mm:ss
#SBATCH --time=0-{args.runtime}:00:00

gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 1
srun gmx mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c final.gro -g production.log -x traj_prod.xtc
"""
    elif args.cluster == "CESGA":
        run_md_script_content = f"""#!/bin/bash -l
#SBATCH -t {args.runtime}:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --gres=gpu:a100

gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 1
srun gmx_mpi mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c confout.gro -g production.log -x traj_prod.xtc
"""
    elif args.cluster == "TETRA": 
        run_md_script_content = f"""#!/bin/bash

#SBATCH --job-name=MD
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -A naiss2023-3-5
#              d-hh:mm:ss
#SBATCH --time=0-{args.runtime}:00:00

gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 1
srun gmx mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c final.gro -g production.log -x traj_prod.xtc"""
    
    run_md_script_path = os.path.join(destination_folder, "run_md.sh")
    with open(run_md_script_path, "w") as run_md_script_file:
        run_md_script_file.write(run_md_script_content)

    # Change permissions to make the script executable
    os.chmod(run_md_script_path, 0o755)

def modify_simulation_time(destination_folder):
    """Modify the simulation time in the prod.mdp file."""
    simulation_time_ns = args.simulation_time
    nsteps = simulation_time_ns * 500000  # 1 ns = 500000 steps

    prod_mdp_path = os.path.join(destination_folder, "prod.mdp")
    with open(prod_mdp_path, "r") as prod_mdp_file:
        lines = prod_mdp_file.readlines()

    with open(prod_mdp_path, "w") as prod_mdp_file:
        for line in lines:
            if line.strip().startswith("nsteps"):
                prod_mdp_file.write(f"nsteps              =  {nsteps}   ; total {simulation_time_ns} ns\n")
            else:
                prod_mdp_file.write(line)

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

# Define the destination folder // TO DO: solve problem when trying to create only FEP from --full_relax false
destination_folder = "md_input_files"

# Call the function to start copying files in the current directory
copy_files_in_directory(".", destination_folder)

# Check if the user wants to prepare FEP files
if args.fep:
    destination_folder_fep = "fep_preparation_files"
    os.makedirs(destination_folder_fep, exist_ok=True)
    prepare_fep_files(args.fep, destination_folder_fep)

    # Define paths to the PDB files
    confout_pdb = os.path.join(destination_folder_fep, "confout_fep.pdb")
    water_pdb = os.path.join(destination_folder_fep, "confout_water.pdb")

    # Process water molecules
    process_water_molecules(confout_pdb, water_pdb)
    clean_pdb(confout_pdb, "L01")
    merge_water_and_pdb(confout_pdb, water_pdb)
    
    # TO DO: Rename and remove elements
    # Rename files with rename_and_remove_elements and remove_last_two_characters
    
    os.remove(water_pdb) # Remove the temporary file of waters

else:
    print("No complex directory specified for FEP preparation.")