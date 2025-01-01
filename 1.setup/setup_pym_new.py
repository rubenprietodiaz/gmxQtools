import os
import argparse
from pymol import cmd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process ligand and protein using PyMOL for PyMemDyn execution.")
    parser.add_argument("--noclean", action='store_true', help="Do not clean the directory after processing.")
    parser.add_argument("-C", "--cluster", choices=["CSB", "CESGA", "TETRA"], default="TETRA", help="Choose the cluster (default: TETRA).")
    parser.add_argument("-p", "--protein", nargs='?', default="protein.pdb", help="Protein file name (default: protein.pdb).")
    parser.add_argument("-l", nargs='?', default="LIG", help="Ligand identifier (default: LIG).")
    parser.add_argument("-cb", "--chargebalance", nargs='?', default="0", help="Charge balance (default: 0).")
    return parser.parse_args()

def clean_protein(protein_file, output_file):
    """Clean the protein file using PyMOL."""
    cmd.load(protein_file, "protein")
    cmd.save(output_file, "protein")
    cmd.delete("all")
    print(f"Protein cleaned and saved as {output_file}")

def create_complex(protein_file, ligand_file, output_file):
    """Combine protein and ligand into a single complex using PyMOL."""
    cmd.load(protein_file, "protein")
    cmd.load(ligand_file, "ligand")
    cmd.save(output_file, "protein or ligand")
    cmd.delete("all")
    print(f"Complex created and saved as {output_file}")

def extract_ligand(complex_file, ligand_identifier, output_file):
    """Extract the ligand from the complex using PyMOL."""
    cmd.load(complex_file, "complex")
    cmd.select("ligand", f"resn {ligand_identifier}")
    cmd.save(output_file, "ligand")
    cmd.delete("all")
    if os.path.getsize(output_file) == 0:
        raise ValueError(f"No ligand atoms found in {complex_file} with identifier {ligand_identifier}")
    print(f"Ligand extracted and saved as {output_file}")

args = parse_arguments()

# Step 1: Clean protein
protein_clean_path = "protein_clean.pdb"
clean_protein(args.protein, protein_clean_path)

# Step 2: Process ligands
for file in os.listdir('.'):
    if file.endswith('.pdb') and file != args.protein and file != protein_clean_path:
        print(f"Processing ligand: {file}")
        dir_name = os.path.splitext(file)[0]
        os.makedirs(dir_name, exist_ok=True)

        # Create complex
        complex_pdb = os.path.join(dir_name, "complex.pdb")
        create_complex(protein_clean_path, file, complex_pdb)

        # Align using PyModSim
        aligned_pdb = os.path.join(dir_name, "complex_aligned.pdb")
        os.system(f'pymodsim -n 3 -p {complex_pdb} > pymodsim.log 2>&1')
        if os.path.exists("finalOutput/homology.pdb"):
            shutil.move("finalOutput/homology.pdb", aligned_pdb)
            print(f"Aligned complex saved as {aligned_pdb}")
        else:
            print("Alignment failed. Check pymodsim.log")
            continue

        # Extract aligned ligand
        aligned_ligand = os.path.join(dir_name, "LIG_aligned.pdb")
        try:
            extract_ligand(aligned_pdb, args.l, aligned_ligand)
        except ValueError as e:
            print(e)
            continue

        # Generate LigParGen parameters
        print(f"Running LigParGen for {dir_name}")
        os.system(f'ligpargen -i {aligned_ligand} -cb {args.chargebalance} -ob 3 -r LIG -n LIG > ligpargen.log 2>&1')

        # Check LigParGen output
        if not os.path.exists("LIG.gmx.gro"):
            print("LigParGen failed, check ligpargen.log")
            continue

        # Rename LigParGen outputs
        shutil.move("LIG.gmx.gro", os.path.join(dir_name, "LIG.gro"))
        shutil.move("LIG.openmm.pdb", os.path.join(dir_name, "LIG.pdb"))
        shutil.move("LIG.gmx.itp", os.path.join(dir_name, "LIG.itp"))
        print(f"LigParGen parameters generated for {aligned_ligand}")

print("All processing complete.")