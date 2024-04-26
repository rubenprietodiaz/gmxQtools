# Execute this bash script in the main directory of the MD simulations
for dir in */ ; do
    cd "$dir"
    echo "Processing $dir"
    mkdir -p pymol
    echo -e "1 0" | gmx trjconv -pbc mol -s topol_prod.tpr -center -ur compact -f traj_prod.xtc -o traj_prod_pymol.xtc
    echo -e "0" | gmx trjconv -s topol_prod.tpr -b 0 -e 0 -o start.gro # First, generate file in gro format
    gmx editconf -f final.gro -o final.pdb
    gmx editconf -f start.gro -o start.pdb # Then convert to pdb to avoid errors
    cp start.pdb pymol/start.pdb
    cp final.pdb pymol/final.pdb
    mv traj_prod_pymol.xtc pymol/traj_prod_pymol.xtc
    echo "All files for $dir are ready for visualization in PyMOL (pymol/ directory)"
    cd ..
done