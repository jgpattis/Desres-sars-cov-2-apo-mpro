#! /usr/bin/env/ python
# Featurize trajectories several different ways to test in next step

import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import os

#traj_num = [ '0000', '0020', '0040', '0060', '0080']
traj_num = [f'{i:04d}' for i in range(100)]
traj_path = '../DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA-'
traj_list = [ traj_path + str(i) + '.dcd' for i in traj_num]
pdb = '../DESRES_protease_chainid.pdb'

# Define features
def filtered_ca_distances(chain=0):
    ''' Pairwize filtered carbon alpha distances defined in filter_distances_01.py'''
    dist_indsA = np.load(open("filtered_dis_ind_10_035_morechainA.npy","rb"))
    dist_indsB = np.load(open("filtered_dis_ind_10_035_morechainB.npy","rb"))
    featurizer = coor.featurizer(pdb)
    if chain == 0:
        featurizer.add_distances(dist_indsA)
    elif chain == 1:
        featurizer.add_distances(dist_indsB)
    else:
        raise ValueError("chain must be 0 or 1")
    return featurizer

def filtered_ca_distances_larger(chain=0):
    ''' Pairwize filtered carbon alpha distances defined in filter_distances_01.py'''
    dist_indsA = np.load(open("filtered_dis_ind_12_03chainA.npy","rb"))
    dist_indsB = np.load(open("filtered_dis_ind_12_03chainA.npy","rb"))
    featurizer = coor.featurizer(pdb)
    if chain == 0:
        featurizer.add_distances(dist_indsA)
    elif chain == 1:
        featurizer.add_distances(dist_indsB)
    else:
        raise ValueError("chain must be 0 or 1")
    return featurizer

def ca_distances_skip5(chain=0):
    ''' Pairwise distances between every 5th carbon alpha '''
    featurizer = coor.featurizer(pdb)
    skip5 = featurizer.select(f'name == CA and chainid == {chain}')[::5]
    featurizer.add_distances(skip5)
    return featurizer

def backbone(chain=0):
    ''' Bachbone Phi and Psi torsion angles '''
    featurizer = coor.featurizer(pdb)
    featurizer.add_backbone_torsions(cossin=True, selstr=f'chainid == {chain}')
    return featurizer

def backbone_chi1(chain=0):
    ''' Bachbone Phi and Psi as well as sidechain chi 1  torsion angles '''
    featurizer = coor.featurizer(pdb)
    featurizer.add_backbone_torsions(cossin=True, selstr=f'chainid == {chain}')
    featurizer.add_sidechain_torsions(which='chi1', cossin=True, selstr=f'chainid == {chain}')
    return featurizer

def backbone_chi1_2(chain=0):
    ''' Bachbone Phi and Psi as well as sidechain chi 1 and chi 2 torsion angles '''
    featurizer = coor.featurizer(pdb)
    featurizer.add_backbone_torsions(cossin=True, selstr=f'chainid == {chain}')
    featurizer.add_sidechain_torsions(which=['chi1','chi2'], cossin=True, selstr=f'chainid == {chain}')
    return featurizer

def sasa_per_res(chain=0):
    ''' Salvent acessable surfase area per residue '''
    def calc_sasa(traj, chain=0, featurizer=None):
        small_traj = traj.atom_slice(atom_indices=featurizer.select(f'chainid == {chain}'))
        res = md.shrake_rupley(small_traj, probe_radius=0.14, n_sphere_points=960, mode='residue')
        return res
    featurizer = coor.featurizer(pdb)
    featurizer.add_custom_func(calc_sasa, dim= int(featurizer.topology.n_residues/2), chain=0, featurizer=featurizer)
    return featurizer

# Add new feature options to feature list and give a 5 character label in feture_short_label
feature_list = (filtered_ca_distances, filtered_ca_distances_larger, ca_distances_skip5, backbone, backbone_chi1, backbone_chi1_2, sasa_per_res)
feature_short_label = ('fil dis', 'm fdis', 'skip5', 'BB Tor', 'BB+C1', 'C1+C2', 'sasa')
chain_list = (0, 1)
label = []
label_file_name = 'feature_list_1.pickl'
stride = 1

# Featurize
for i,j in enumerate(feature_list):
    for k in chain_list:
        file_name = f"feature_data_02/{j.__name__}_chain_{k}.h5"
        feat = j(chain=k)
        if k == 1:
            label.append({
                'long_label': j.__name__,
                'short_label': feature_short_label[i],
                'path_a': f'feature_data_02/{j.__name__}_chain_0.h5',
                'path_b': f'feature_data_02/{j.__name__}_chain_1.h5',
                'num_features': feat.dimension()
            })
        if os.path.exists(file_name):
            print(f'{j.__name__} exists')
            continue
        reader = coor.source(traj_list, features=feat, stride=stride)
        reader.write_to_hdf5(file_name)

# Save as pandas dataframe
import pandas as pd
results = pd.DataFrame(label)
print(results.head())
results.to_pickle(label_file_name)


