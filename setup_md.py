import os
import shutil
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to copy files, create submit.sh, and modify simulation duration.")
    parser.add_argument("-d", "--output-dir", default="md_input_files", help="Output directory (default: md_input_files).")
    parser.add_argument("-t", "--simulation-time", type=int, help="Simulation time in nanoseconds.")
    # Check README to see what I have to add here
    return parser.parse_args()

def copy_files_in_directory(directory, destination_folder):
    # Iterate over the folders within the current directory
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            current_folder = os.path.join(root, dir)
            # Call the copy_files function to copy files in the current folder
            copy_files(current_folder, destination_folder)

def copy_files(folder, destination_folder):
    # Check the existence of each required file
    if all(os.path.isfile(os.path.join(folder, file)) for file in ["prod.mdp", "topol.top", "index.ndx"]):
        # Extract the name of the current folder
        folder_name = os.path.basename(folder)

        # Create the folder for files within the destination folder
        destination_folder_path = os.path.join(destination_folder, folder_name)
        os.makedirs(destination_folder_path, exist_ok=True)

        # Copy the files to the corresponding folder
        for file in ["prod.mdp", "topol.top", "index.ndx"]:
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

        # Create the submit.sh script inside the destination folder
        create_submit_script(destination_folder_path)

        # Create the run_md.sh script inside the destination folder
        create_run_md_script(destination_folder_path)

        # Modify prod.mdp if simulation time is provided
        if args.simulation_time:
            modify_simulation_time(destination_folder_path)

def create_submit_script(destination_folder):
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
    run_md_script_content = """#!/bin/bash

#SBATCH --job-name=MD
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -A naiss2023-3-5
#              d-hh:mm:ss
#SBATCH --time=0-24:00:00

gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 1
srun gmx mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c final.gro -g production.log -x traj_prod.xtc
"""
    run_md_script_path = os.path.join(destination_folder, "run_md.sh")
    with open(run_md_script_path, "w") as run_md_script_file:
        run_md_script_file.write(run_md_script_content)

    # Change permissions to make the script executable
    os.chmod(run_md_script_path, 0o755)

def modify_simulation_time(destination_folder):
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

# Parse arguments
args = parse_arguments()

# Define the destination folder
destination_folder = args.output_dir

# Call the function to start copying files in the current directory
copy_files_in_directory(".", destination_folder)

# Create the submit.sh script inside the destination folder
create_submit_script(destination_folder)

print(f"Folder {destination_folder} with scripts and necessary files created. Run 'cd {destination_folder}' and 'bash submit.sh' to start the simulations.")
