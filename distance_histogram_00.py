#! /usr/bin/env python
# make histogram to determine CA distances with larger minimum distance to be filtered out
# make histogram to determine CA distances with smaller standard deviations to be filtered out

import mdtraj as md
import pyemma.coordinates as coor
import matplotlib.pyplot as plt
import numpy as np

feat = coor.featurizer('../DESRES_protease_chainid.pdb')
traj_num = [f'{i:04d}' for i in range(100)] 
traj_path = '../DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA-'
traj_list = [ traj_path + str(i) + '.dcd' for i in traj_num]

feat.add_distances(feat.pairs(feat.select('name == CA and chainid == 0'), excluded_neighbors=3))
traj = coor.load(traj_list, feat, stride=5)

traj_cat = np.concatenate(traj)
min_dist = traj_cat.min(axis=0)
std_dist = traj_cat.std(axis=0)

fig1, ax1 = plt.subplots()
ax1.hist(min_dist, bins=60)
ax1.set_xlabel('Minimum Distance', fontsize=16)
ax1.set_ylabel('Frequency', fontsize=16)
fig1.savefig('histogram_00/distance_histogram_1.pdf')

fig2, ax2 = plt.subplots()
ax2.hist(std_dist, bins=60)
ax2.set_xlabel('Distance Standard Deviation', fontsize=16)
ax2.set_ylabel('Frequency', fontsize=16)
fig2.savefig('histogram_00/standard_deviation_histogram_1.pdf')

