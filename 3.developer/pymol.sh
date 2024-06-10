# This script is used to generate the necessary files for visualization in PyMOL
# Execute this bash script in the main directory of the MD simulations
for dir in */ ; do
    cd "$dir"
    echo "Processing $dir"
    mkdir -p pymol
    echo -e "1 0" | gmx trjconv -pbc mol -s topol_prod.tpr -center -ur compact -f traj_prod.xtc -o traj_prod_pymol.xtc &>> visualization.log
    echo -e "0" | gmx trjconv -s topol_prod.tpr -b 0 -e 0 -o start.gro &>> visualization.log # First, generate file in gro format
    gmx editconf -f final.gro -o final.pdb &>> visualization.log
    gmx editconf -f start.gro -o start.pdb  &>> visualization.log # Then convert to pdb to avoid errors
    cp start.pdb pymol/start.pdb
    cp final.pdb pymol/final.pdb
    mv traj_prod_pymol.xtc pymol/traj_prod_pymol.xtc

    # Create the PyMOL script // TO-DO: IMPROVE THIS (some problems with selection and ini state) // Use template file?
    cat > pymol/load_md.pml << EOF
#!/usr/bin/env python
set defer_builds_mode, 3
load start.pdb, ini-state
load start.pdb, equi
color grey70, ini-state
hide everything, resn POPC
hide lines, (element h and neighbor element c)
select protein, chain A and equi
select solvent, resn SOL
select membrane, resn POPC
select memblimi, equi and name n4+p8+na*
hide everything, membrane
hide everything, solvent
show cartoon, protein
show spheres, solvent and elem O
show spheres, memblimi
set sphere_scale, 0.2
cmd.spectrum("count",selection="(protein)&e. c")

select ligand, equi and resn lig
show sticks, ligand and not (hydro and neighbor element c)
util.cba(154,"ligand",_self=cmd)
cmd.disable('ligand')

load_traj traj_prod_pymol.xtc, equi
intra_fit equi and name ca

center ligand
zoom ligand
set auto_zoom, 0
deselect
EOF

    echo "All files are ready for visualization in PyMOL (execute pymol load_md.pml)"
    cd ..
done