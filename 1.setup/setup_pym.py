import os
import argparse
import shutil

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process ligand and protein for PyMemDyn execution.")
    parser.add_argument("--noclean", action='store_true', help="Do not clean the directory after processing.")
    parser.add_argument("-C", "--cluster", choices=["CSB", "CESGA", "TETRA"], default="TETRA", help="Choose the cluster (default: TETRA).") # Changed to CSB for testing
    parser.add_argument("-p", "--protein", nargs='?', default="protein.pdb", help="Protein file name (default: protein.pdb).")
    parser.add_argument("-l", nargs='?', default="LIG", help="Ligand identifier (default: LIG).")
    parser.add_argument("-na", "--noalign", action='store_true', help="Do not align the protein and ligand. Use this option if the complex is already aligned with PyModSim.")

    # Parameters for LigParGen execution
    parser.add_argument("-cb", "--chargebalance", nargs='?', default="0", help="Charge balance (default: 0).") # Test with 0 and 1 (Lundbeck ligands with piperazine ring)

    # Parameters for Pymemdyn execution
    parser.add_argument("-r", "--res", nargs='?', default="ca", help="Restraints. Options: bw (Ballesteros-Weinstein Restrained Relaxation), ca (C-Alpha Restrained Relaxation - default).")
    parser.add_argument("-w", nargs='?', default="HOH", help="Water identifiers (default: HOH).")
    parser.add_argument("-i", nargs='?', default="NA", help="Ion identifiers (default: NA).")
    parser.add_argument("-fep", "--fep", action='store_true', help="Choose to prepare files for FEP calculations (add --full_relax false to pymemdyn script).")
    return parser.parse_args()

args = parse_arguments()
start_dir = os.getcwd()
if args.fep:
    print('You have chosen to prepare files for FEP calculations. This is faster, but you will need to run the full relaxation manually if want to run MD.')
if args.noclean:
    print('You have chosen to not clean the directory after processing')
if args.cluster:
    print(f'You have chosen the {args.cluster} cluster for execution. Please, make sure the options are correct.')
if args.noalign:
    print(f'You have chosen to avoid PyModSim alignment. Be sure the complex is already aligned.')
# Loop through all .pdb files in the current directory, excluding the protein file
for file in os.listdir('.'):
    if file.endswith('.pdb') and file != args.protein:
        print("Processing:", file)

        # Read the content of the file, remove unnecessary lines and replace UNK with LIG
        with open(file, 'r') as f:
            content = f.readlines()
        content = [line for line in content if not line.startswith(('END', 'CONECT', 'REMARK', 'TITLE'))]
        content = [line.replace(args.l, 'LIG') for line in content]

        # Create a directory named after the ligand file (without the .pdb extension) and write the modified ligand file named as LIG.pdb
        dir_name = os.path.splitext(file)[0]
        os.makedirs(dir_name, exist_ok=True)
        ligand_path = os.path.join(dir_name, 'LIG.pdb')
        with open(ligand_path, 'w') as f:
            f.writelines(content)

        # Concatenate protein file and ligand file to create complex.pdb
        protein_path = args.protein
        complex_path = os.path.join(dir_name, 'complex.pdb')
        with open(protein_path, 'r') as f_protein, open(ligand_path, 'r') as f_ligand, open(complex_path, 'w') as f_complex:
            f_complex.write(f_protein.read())
            f_complex.write(f_ligand.read())

        # Enter the directory
        os.chdir(dir_name)

        # Execute ligpargen and rename files
        print(f"[1/3] Running LigParGen for {dir_name}")
        os.system(f'ligpargen -i LIG.pdb -cb {args.chargebalance} -ob 3 -r LIG -n LIG > ligpargen.log 2>&1')
        os.rename('LIG.gmx.gro', 'LIG.gro')
        os.rename('LIG.openmm.pdb', 'LIG.pdb')
        os.rename('LIG.gmx.itp', 'LIG.itp')

        if args.noalign is False:
            # Execute pymodsim for alignment of complex.pdb and clean files
            print(f"[2/3] Running PyModSim for {dir_name}")
            os.system('pymodsim -n 3 -p complex.pdb > pymodsim.log 2>&1')
            if os.path.exists('finalOutput/complex.pdb'):
                os.rename('complex.pdb', 'complex.pdb.bak') # Backup the original complex.pdb
                shutil.copy('finalOutput/complex.pdb', 'complex.pdb')
                
                # Remove files from pymodsim
                if args.noclean is False:
                    files_to_delete = [
                        'pymodsim.log',
                        'Model_output.tgz',
                        'protein_stripped.pdb',
                        '2membranes.inp',
                        'datapar2',
                        'protein_strippedout.pdb',
                        'fort.41',
                        'fort.4',
                        'datapar1',
                        'datasub1',
                        'protein_aligned.pdb',
                        'homology.pdb'
                    ]
                    for file in files_to_delete:
                        if os.path.exists(file):
                            os.remove(file)
                    os.system('rm -rf finalOutput/')       
            else:
                print("Warning: PyModSim alignment didn't work, check manually the complex.pdb file.")
        else:
            print(f"[2/3] Skipping alignment for {dir_name}.")
        # pymemdyn.sh script inside the directory
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
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=ruben.prieto@usc.es
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
        print(f"[3/3] pymemdyn.sh script ({args.cluster}) created for {dir_name}")
        os.chdir(start_dir)

        print(f"Processing of {dir_name} complete.")

# Backup the provided PDB files
os.makedirs('inputFiles', exist_ok=True)
os.system('mv *.pdb inputFiles/')

# submit_pym.sh script
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
print(f'All ligands processed. Provided PDB files were moved to inputFiles directory.')
print(f'Execute submit_pym.sh in {args.cluster} to start pymemdyn equilibration.')
print(f'After running, you can use setup_md.py or setup_fep.py to prepare the files for MD or FEP calculations respectively.')
print(f'Check gmxQtools repository for further information.')