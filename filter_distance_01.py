#! /usr/bin/env/ python
# filter out CA distances with large minimum
# filter out CA distances with small standard deviations
# save to file for later use
# will plot the distances used

import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
from plot_structure_util import plot_vmd_cylinder_from_inds, plot_pymol_cylinder_from_inds

dis_cutoff = 1.0
std_cutoff = 0.035
outfile = 'filtered_dis_ind_10_035_more'
save = True
plot = 'all'  # should be all, pymol, vmd, or none

traj_num = [f'{i:04d}' for i in range(100)] 
traj_path = '../DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA-'
traj_list = [ traj_path + str(i) + '.dcd' for i in traj_num]

pdb = '../DESRES_protease_chainid.pdb'
feat = coor.featurizer(pdb)
feat.add_distances(feat.pairs(feat.select('name == CA and chainid == 0'), excluded_neighbors=3))
traj = coor.load(traj_list, feat, stride=5)
traj_cat = np.concatenate(traj)

feat1 = coor.featurizer(pdb)
feat1.add_distances(feat1.pairs(feat1.select('name == CA and chainid == 1'), excluded_neighbors=3))
traj1 = coor.load(traj_list, feat, stride=5)
traj_cat1 = np.concatenate(traj)

traj_cat_pair = np.concatenate((traj_cat, traj_cat1), axis=0)

min_dist = traj_cat_pair.min(axis=0)
std_dist = traj_cat_pair.std(axis=0)

new_dists = np.where((min_dist < dis_cutoff) & (std_dist > std_cutoff))[0]

print('new distances:', new_dists.shape)

if save == True:
    out = np.zeros((len(new_dists), 2), dtype=np.int16)
    label = feat.describe()
    for i,j in enumerate(new_dists):
        tmp = label[j].split()
        out[i,0] = int(tmp[4])
        out[i,1] = int(tmp[10])

    with open(outfile+'chainA.npy','wb') as handle:
        np.save(handle, out)

    out1 = np.zeros((len(new_dists), 2), dtype=np.int16)
    label1 = feat1.describe()
    for i,j in enumerate(new_dists):
        tmp1 = label1[j].split()
        out1[i,0] = int(tmp1[4])
        out1[i,1] = int(tmp1[10])

    with open(outfile+'chainB.npy','wb') as handle:
        np.save(handle, out1)

    if plot == 'all':
        plot_vmd_cylinder_from_inds(pdb, out, outfile + '_draw_cylinder')
        print(plot_vmd_cylinder_from_inds.__doc__)
        plot_pymol_cylinder_from_inds(pdb, out, outfile + '_draw_cylinder')
        print(plot_pymol_cylinder_from_inds.__doc__)
    elif plot == 'vmd':
        plot_vmd_cylinder_from_inds(pdb, out, outfile + '_draw_cylinder')
        print(plot_vmd_cylinder_from_inds.__doc__)
    elif plot == 'pymol':
        plot_pymol_cylinder_from_inds(pdb, out, outfile + '_draw_cylinder')
        print(plot_pymol_cylinder_from_inds.__doc__)
