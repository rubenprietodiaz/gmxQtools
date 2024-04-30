import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process ligand and protein for PyMemDyn execution.")
    parser.add_argument("-C", "--cluster", choices=["CSB", "CESGA", "TETRA"], default="TETRA", help="Choose the cluster (default: TETRA).")
    parser.add_argument("-p", "--protein", nargs='?', default="protein.pdb", help="Protein file name (default: protein.pdb).")
    parser.add_argument("-l", nargs='?', default="LIG", help="Ligand identifier (default: LIG).")

    # Parameters for Pymemdyn execution
    parser.add_argument("-r", "--res", nargs='?', default="ca", help="Restraints. Options: bw (Ballesteros-Weinstein Restrained Relaxation), ca (C-Alpha Restrained Relaxation - default).")
    parser.add_argument("-w", nargs='?', default="HOH", help="Water identifiers (default: HOH).")
    parser.add_argument("-i", nargs='?', default="NA", help="Ion identifiers (default: NA).")
    parser.add_argument("-fep", "--fep", action='store_true', help="Choose to prepare files for FEP calculations (add --full_relax false to pymemdyn script).")
    return parser.parse_args()

# Parse arguments
args = parse_arguments()

# Save the current directory
start_dir = os.getcwd()

# Loop through all .pdb files in the current directory, excluding the protein file
for file in os.listdir('.'):
    if file.endswith('.pdb') and file != args.protein:
        print("Processing:", file)
        
        # Read the content of the file
        with open(file, 'r') as f:
            content = f.readlines()
        
        # Remove unnecessary lines and replace UNK with LIG
        content = [line for line in content if not line.startswith(('END', 'CONECT', 'REMARK', 'TITLE'))]
        content = [line.replace(args.l, 'LIG') for line in content]
        
        # Create a directory named after the ligand file (without the .pdb extension)
        dir_name = os.path.splitext(file)[0]
        os.makedirs(dir_name, exist_ok=True)
        
        # Write the modified ligand file to the new directory and rename it to LIG.pdb
        ligand_path = os.path.join(dir_name, 'LIG.pdb')
        with open(ligand_path, 'w') as f:
            f.writelines(content)
        
        # Concatenate protein file and ligand file to create complex.pdb
        protein_path = args.protein
        complex_path = os.path.join(dir_name, 'complex.pdb')
        with open(protein_path, 'r') as f_protein, open(ligand_path, 'r') as f_ligand, open(complex_path, 'w') as f_complex:
            f_complex.write(f_protein.read())
            f_complex.write(f_ligand.read())
        
        # Change into the ligand's directory
        os.chdir(dir_name)
        
        # Execute ligpargen
        os.system('ligpargen -i LIG.pdb -cb 0 -ob 3 -r LIG -n LIG')
        
        # Rename files
        os.rename('LIG.gmx.gro', 'LIG.gro')
        os.rename('LIG.openmm.pdb', 'LIG.pdb')
        os.rename('LIG.gmx.itp', 'LIG.itp')
        
        # Create pymemdyn.sh script
        if args.cluster == "CSB":
            pymemdyn_content = f"""#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 24:00:00
#SBATCH --gpus-per-task=1
#SBATCH --job-name=pymemdyn
pymemdyn -p complex.pdb --res {args.res} -w {args.w} -i {args.i} -l LIG {'--full_relax false' if args.fep else ''}\n"""
        elif args.cluster == "CESGA":
            pymemdyn_content = f"""#!/bin/bash -l
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem-per-cpu=4G
#SBATCH -t 24:00:00
#SBATCH --job-name=pymemdyn
pymemdyn -p complex.pdb --res {args.res} -w {args.w} -i {args.i} -l LIG {'--full_relax false' if args.fep else ''}\n"""
        elif args.cluster == "TETRA":
            pymemdyn_content = f"""#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 24:00:00
#SBATCH --job-name=pymemdyn
pymemdyn -p complex.pdb --res {args.res} -w {args.w} -i {args.i} -l LIG {'--full_relax false' if args.fep else ''}\n"""

        
        with open('pymemdyn.sh', 'w') as f_pymemdyn:
            f_pymemdyn.write(pymemdyn_content)
        
        # Return to the start directory
        os.chdir(start_dir)
        
        print("Created complex.pdb, executed ligpargen, renamed files, and created pymemdyn.sh in", dir_name)

# Create directory for input files and move all *.pdb to that directory
os.makedirs('inputFiles', exist_ok=True)
os.system('mv *.pdb inputFiles/')
print("Created backup of input files in inputFiles/.")

# Create submit_pym.sh to submit jobs in all directories
with open('submit_pym.sh', 'w') as f_submit:
    f_submit.write("#!/bin/bash\n\n")
    f_submit.write("echo 'Processing directories:'\n")
    f_submit.write("start_dir=$(pwd)\n")
    f_submit.write("for folder in ./*; do\n")
    f_submit.write("    if [ -d \"$folder\" ]; then\n")
    f_submit.write("        cd \"$folder\" || continue\n")
    f_submit.write("        sbatch pymemdyn.sh\n")
    f_submit.write("        cd \"$start_dir\"\n")
    f_submit.write("    fi\n")
    f_submit.write("done\n")
    f_submit.write("echo 'All jobs submitted.'\n")

print("submit_pym.sh created successfully.")
print("All processes complete. Pymemdyn execution ready, run 'sh submit_pym.sh' to start.")
