import os
import shutil
import argparse
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description="Setup your md simulation after pymemdyn equilibration.")
    parser.add_argument("-t", "--simulation-time", type=int, help="Simulation time in nanoseconds (default by PyMemDyn: 10 ns).")
    parser.add_argument("-rt", "--runtime", type=int, help="Runtime in hours (default: 24).", default=24)
    parser.add_argument("-C", "--cluster", choices=["CSB", "CESGA", "TETRA"], default="TETRA", help="Choose the cluster (default: TETRA).")
    parser.add_argument("-n", "--num-replicas", type=int, help="Number of replicas for the MD simulations (default: 1).", default=1)
    return parser.parse_args()

def copy_files_in_directory(directory, destination_folder):
    """Iterate over the folders within the current directory and call the copy_files function."""
    # Iterate over the folders within the current directory
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            current_folder = os.path.join(root, dir)
            # Call the copy_files function to copy files in the current folder
            for replica in range(args.num_replicas):
                replica_folder = os.path.join(destination_folder, f"{dir}_{replica+1}")
                copy_files(current_folder, replica_folder)
    # Create submit script outside the loop
    create_submit_script(destination_folder)
    print(f"Folder {destination_folder} with scripts and necessary files created. Run 'cd {destination_folder}' and 'sh submit_md.sh' to start the simulations.")

def copy_files(folder, destination_folder):
    """Copy the required files for MD simulations to the destination folder."""
    # Check the existence of each required file
    if all(os.path.isfile(os.path.join(folder, file)) for file in ["prod.mdp", "topol.top", "index.ndx", "topol.tpr"]) and os.path.exists(os.path.join(folder, "finalOutput", "confout.gro")):
        print(f"All the required files are present in the folder {folder}.")
        # Create the folder for files within the destination folder
        os.makedirs(destination_folder, exist_ok=True)

        # Copy the files to the corresponding folder
        for file in ["prod.mdp", "topol.top", "index.ndx", "topol.tpr"]:
            shutil.copy(os.path.join(folder, file), destination_folder)

        # Check the existence of the finalOutput folder
        final_output_folder = os.path.join(folder, "finalOutput")
        if os.path.isdir(final_output_folder):
            confout_gro = os.path.join(final_output_folder, "confout.gro")
            if os.path.isfile(confout_gro):
                shutil.copy(confout_gro, destination_folder)

        # Copy .itp files
        for itp_file in os.listdir(folder):
            if itp_file.endswith(".itp"):
                shutil.copy(os.path.join(folder, itp_file), destination_folder)
        
        print(f"Files copied to {destination_folder}.")

        # Create the run_md.sh script inside each destination folder
        create_run_md_script(destination_folder)

        # Modify prod.mdp if simulation time is provided, if not, prod.mdp will be used as is (10 ns)
        if args.simulation_time:
            modify_simulation_time(destination_folder)
        
        # Modify gen_seed for each replica
        modify_gen_seed(destination_folder)
    else:
        return False

def create_submit_script(destination_folder):
    """Create the submit_md.sh script to submit the MD simulations."""
    submit_script_content = """#!/bin/bash

start_dir=$(pwd)

for dir in */ ; do
    cd "$dir"
    sbatch run_md.sh
    cd "$start_dir/$1"
done
"""
    submit_script_path = os.path.join(destination_folder, "submit_md.sh")
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
#SBATCH --time=0-{args.runtime}:00:00

module load gromacs

gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 1
srun gmx mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c final.gro -g production.log -x traj_prod.xtc

# Post-processing for visualization
mkdir -p finalOutput
echo -e "1 0" | gmx trjconv -pbc mol -s topol_prod.tpr -center -ur compact -f traj_prod.xtc -o traj_prod_pymol.xtc &>> visualization.log
echo -e "0" | gmx trjconv -s topol_prod.tpr -b 0 -e 0 -o start.gro &>> visualization.log
gmx editconf -f final.gro -o final.pdb &>> visualization.log
gmx editconf -f start.gro -o start.pdb &>> visualization.log

cp start.pdb finalOutput/start.pdb
cp final.pdb finalOutput/final.pdb
mv traj_prod_pymol.xtc finalOutput/traj_prod_pymol.xtc
"""
    elif args.cluster == "CESGA":
        run_md_script_content = f"""#!/bin/bash -l
#SBATCH -t {args.runtime}:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --gres=gpu:a100

module load gromacs

gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 1
srun gmx_mpi mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c confout.gro -g production.log -x traj_prod.xtc

# Post-processing for visualization
mkdir -p finalOutput
echo -e "1 0" | gmx trjconv -pbc mol -s topol_prod.tpr -center -ur compact -f traj_prod.xtc -o traj_prod_pymol.xtc &>> visualization.log
echo -e "0" | gmx trjconv -s topol_prod.tpr -b 0 -e 0 -o start.gro &>> visualization.log
gmx editconf -f final.gro -o final.pdb &>> visualization.log
gmx editconf -f start.gro -o start.pdb &>> visualization.log

cp start.pdb finalOutput/start.pdb
cp final.pdb finalOutput/final.pdb
mv traj_prod_pymol.xtc finalOutput/traj_prod_pymol.xtc
"""
    elif args.cluster == "TETRA":
        run_md_script_content = f"""#!/bin/bash

#SBATCH --job-name=MD
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -A naiss2023-3-5
#SBATCH --time=0-{args.runtime}:00:00

module load gromacs

gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 1
srun gmx mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c final.gro -g production.log -x traj_prod.xtc

# Post-processing for visualization
mkdir -p finalOutput
echo -e "1 0" | gmx trjconv -pbc mol -s topol_prod.tpr -center -ur compact -f traj_prod.xtc -o traj_prod_pymol.xtc &>> visualization.log
echo -e "0" | gmx trjconv -s topol_prod.tpr -b 0 -e 0 -o start.gro &>> visualization.log
gmx editconf -f final.gro -o final.pdb &>> visualization.log
gmx editconf -f start.gro -o start.pdb &>> visualization.log

cp start.pdb finalOutput/start.pdb
cp final.pdb finalOutput/final.pdb
mv traj_prod_pymol.xtc finalOutput/traj_prod_pymol.xtc
"""
    run_md_script_path = os.path.join(destination_folder, "run_md.sh")
    with open(run_md_script_path, "w") as run_md_script_file:
        run_md_script_file.write(run_md_script_content)

    # Make the script executable
    os.chmod(run_md_script_path, 0o755)

def modify_simulation_time(destination_folder):
    """Modify the simulation time in the prod.mdp file."""
    simulation_time_ns = args.simulation_time
    nsteps = simulation_time_ns * 500000  # 1 ns = 500000 steps

    prod_mdp_path = os.path.join(destination_folder, "prod.mdp")
    with open(prod_mdp_path, "r") as prod_mdp_file:
        lines = prod_mdp_file.readlines()

    with open(prod_mdp_path, "w") as prod_mdp_file: # Add a function to change gen_seed (if the runs start at exactly the same time, they will have the same seed)
        for line in lines:
            if line.strip().startswith("nsteps"):
                prod_mdp_file.write(f"nsteps              =  {nsteps}   ; total {simulation_time_ns} ns\n")
            else:
                prod_mdp_file.write(line)

def modify_gen_seed(destination_folder):
    """Modify the gen_seed in the prod.mdp file to ensure different random seeds for each replica."""
    import random
    gen_seed = random.randint(1, 2147483647)  # Random seed between 1 and 2147483647

    prod_mdp_path = os.path.join(destination_folder, "prod.mdp")
    with open(prod_mdp_path, "r") as prod_mdp_file:
        lines = prod_mdp_file.readlines()

    with open(prod_mdp_path, "w") as prod_mdp_file:
        for line in lines:
            if line.strip().startswith(";gen_seed"):
                prod_mdp_file.write(f"gen_seed            =  {gen_seed}\n")
            else:
                prod_mdp_file.write(line)
        # # If gen_seed is not present, add it at the end
        # if not any(line.strip().startswith(";gen_seed") for line in lines):
        #     prod_mdp_file.write(f"gen_seed            =  {gen_seed}\n")

# RUN
args = parse_arguments()
destination_folder = "md_input_files"
copy_files_in_directory(".", destination_folder)
