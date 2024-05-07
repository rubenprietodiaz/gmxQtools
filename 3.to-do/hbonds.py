import mdtraj as md

# Define input files
pdb_file = 'start.pdb'
traj_file = 'traj_prod.xtc'

# Define selection criteria
selection = 'resn L01 or protein'

traj = md.load(traj_file, top=pdb_file)
print('Trajectory loaded successfully.')
traj_slice = traj.topology.select(selection)
print('Trajectory sliced successfully.')
subtraj = traj.atom_slice(traj_slice)
print('Subtrajectory sliced successfully.')

hbonds = md.baker_hubbard(subtraj, freq=0.0)

for hbond in hbonds:
    donor, hydrogen, acceptor = hbond[:]

    donor_name = str(subtraj.topology.atom(donor))[7:]
    acceptor_name = str(subtraj.topology.atom(acceptor))[7:]

    donor_res_index = subtraj.topology.atom(donor).residue.index
    donor_res = subtraj.topology.residue(donor_res_index)
    donor_res_name = donor_res.name

    acceptor_res_index = subtraj.topology.atom(acceptor).residue.index
    acceptor_res = subtraj.topology.residue(acceptor_res_index)
    acceptor_res_name = acceptor_res.name

    if donor_res_name == 'L01' or acceptor_res_name == 'L01':
        print(donor_name, donor_res_name, donor_res_index, acceptor_name, acceptor_res_name, acceptor_res_index)
