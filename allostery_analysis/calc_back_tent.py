import mdentropy.metrics
import mdtraj as md
import numpy as np
import pickle

outname = 'dtent'
first_residue = 1 # residue number of the first residue in the protein

#traj_num = [ '0000', '0020', '0040', '0060', '0080']
traj_num = [f'{i:04d}' for i in range(100)]
traj_path = '../../DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA-'
traj_list = [ traj_path + str(i) + '.dcd' for i in traj_num]
pdb = '../../DESRES_protease_chainid.pdb'

traj = md.load(traj_list, stride=20, top=pdb)
traj1 = traj[:-1]
traj2 = traj[1:]

dte = mdentropy.metrics.DihedralTransferEntropy(types=['phi', 'psi'], n_bins=3, method='knn', normed=False)
out = dte.partial_transform((traj1, traj2), shuffle=5)

print('entropy calulation is complete for backbone dihedral')

np.save(f'output/{outname}.npy', out)

np.save(f'output/protease_residue_labels_{outname}.npy', dte.labels + first_residue)
