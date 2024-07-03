# Install this script in .pymol/scripts folder as cyx.pml
# Add 'alias cyx, run ~/.pymol/scripts/cyx.pml'
# Before QligFEP setup.py, execute cyx in pymol with protein.pdb (after protprep.py in QligFEP protocol)
# You will see the SG atom numbers for CYX bonds

select cyx, resn CYX
select sulphur, cyx and elem S
show lines, sulphur
hide everything, not sulphur
label sulphur, ID