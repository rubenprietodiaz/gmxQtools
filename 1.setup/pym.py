import os
import argparse
import shutil
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process ligand and protein for PyMemDyn execution.")
    parser.add_argument("--noclean", action='store_true', help="Do not clean the directory after processing.")
    parser.add_argument("-C", "--cluster", choices=["CSB", "CESGA", "TETRA"], default="TETRA", help="Choose the cluster (default: TETRA).")
    parser.add_argument("-p", "--protein", nargs='?', default="protein.pdb", help="Protein file name (default: protein.pdb).")
    parser.add_argument("-l", nargs='?', default="LIG", help="Ligand identifier (default: LIG).")
    parser.add_argument("-na", "--noalign", action='store_true', help="Do not align the protein and ligand. Use this option if the complex is already aligned with PyModSim.")
    parser.add_argument("-cb", "--chargebalance", nargs='?', default="0", help="Charge balance (default: 0).")
    parser.add_argument("-r", "--res", nargs='?', default="ca", help="Restraints. Options: bw (Ballesteros-Weinstein Restrained Relaxation), ca (C-Alpha Restrained Relaxation - default).")
    parser.add_argument("-w", nargs='?', default="HOH", help="Water identifiers (default: HOH).")
    parser.add_argument("-i", nargs='?', default="NA", help="Ion identifiers (default: NA).")
    parser.add_argument("-fep", "--fep", action='store_true', help="Choose to prepare files for FEP calculations (add --full_relax false to pymemdyn script).")
    return parser.parse_args()

args = parse_arguments()
start_dir = os.getcwd()
if args.fep:
    print('You have chosen to prepare files for FEP calculations. This is faster, but you will need to run the full relaxation manually if you want to run MD.')
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

        # Create a directory named after the ligand file (without the .pdb extension)
        dir_name = os.path.splitext(file)[0]
        os.makedirs(dir_name, exist_ok=True)

        # Paths for input/output files
        protein_path = args.protein
        ligand_path = os.path.join(dir_name, 'LIG.pdb')
        complex_path = os.path.join(dir_name, 'complex.pdb')

        # Copy the ligand file into the directory
        shutil.copy(file, ligand_path)

        # Generate complex.pdb using PyMOL (CMD mode)
        print(f"[1/4] Generating complex.pdb for {dir_name} using PyMOL CMD mode")
        pymol_cmd = f"""
from pymol import cmd
cmd.load("{protein_path}", "protein")
cmd.load("{ligand_path}", "ligand")
cmd.alter("ligand", "resn='{args.l}'")
cmd.save("{complex_path}", "protein or ligand")
cmd.quit()
"""
        try:
            with open("create_complex.pml", "w") as script_file:
                script_file.write(pymol_cmd)
            subprocess.run(["pymol", "-cq", "create_complex.pml"], check=True)
            os.remove("create_complex.pml")
        except subprocess.CalledProcessError:
            print(f"Error: PyMOL CMD failed to generate complex.pdb for {dir_name}. Check your input files.")
            continue

        os.remove(ligand_path)  # Remove the ligand file from the directory before alignment to avoid confusion

#         Extract the ligand from complex.pdb using PyMOL (for testing purposes)
#         print(f"[2/4] Extracting ligand from complex.pdb for {dir_name}")
#         pymol_extract_cmd = f"""
# from pymol import cmd
# cmd.load("{complex_path}", "complex")
# cmd.select("ligand", "resn {args.l}")
# cmd.save("{ligand_path}", "ligand")
# cmd.quit()
# """
#         try:
#             with open("extract_ligand.pml", "w") as extract_script:
#                 extract_script.write(pymol_extract_cmd)
#             subprocess.run(["pymol", "-cq", "extract_ligand.pml"], check=True)
#             os.remove("extract_ligand.pml")
#             print(f"Ligand extracted successfully: {ligand_path}")
#         except subprocess.CalledProcessError:
#             print(f"Error: Failed to extract ligand from complex.pdb for {dir_name}. Check the complex.pdb structure.")

        # Enter the directory
        os.chdir(dir_name)
        
        # Execute pymodsim for alignment of complex.pdb and clean files
        print(f"[2/4] Running PyModSim for {dir_name}")
        os.system('pymodsim -n 3 -p complex.pdb > pymodsim.log 2>&1')
        if os.path.exists('finalOutput/homology.pdb'):
            os.rename('complex.pdb', 'complex.pdb.bak')  # Backup the original complex.pdb
            shutil.copy('finalOutput/homology.pdb', 'complex.pdb')    
        else:
            print("Warning: PyModSim alignment didn't work, check manually the complex.pdb file.")

