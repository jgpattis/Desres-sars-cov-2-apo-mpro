#! /usr/bin/env/ python
# Featurize trajectories several different ways to test in next step

import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import os
import pandas as pd

var_cutoff=0.9 # adjust to find elbow of of cumulative kinetic variance

feature_list = pd.read_pickle('feature_list_1.pickl')

data_a = coor.load('feature_data_02/backbone_chi1_2_chain_0.h5')
data_b = coor.load('feature_data_02/backbone_chi1_2_chain_1.h5')
data = data_a + data_b

tica = coor.tica(data=data, lag=10, kinetic_map=False, commute_map=True)

with open('tica_data_05/chi2_cumvar.npy','wb') as handle:
    np.save(handle, tica.cumvar)

tica.var_cutoff = var_cutoff

print('Number of dimentions saved is: ', tica.dimension())
tica.write_to_hdf5('tica_data_05/chi2_tica_data.h5')
