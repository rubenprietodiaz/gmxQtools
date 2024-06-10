# TO-DO
- Add function to change random seed in setup_md.py
- Add rsync function after setup_pym.py to send jobs to cluster. Indicate the directory with an argument.
- Include pymol.sh in md_setup and execute after the end of MD simulation (then, all the compounds will have pymol directory inside)

# Solved
- Correct problem with chains in setup_fep.py (change HOH A to HOH X may be a solution): SOLVED, renaming water chains to W
- Add pymodsim to setup_pym.py (for each complex.pdb, after execution of pymodsim, take finalOutput/complex.pdb) -> Only in dev