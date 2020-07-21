import pyemma.coordinates as coor
import numpy as np
import pyemma.msm as msm
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt

sys = 'fdis'
n_clusters = 100
dtrajs = coor.load(f'cluster_data_10/{sys}_{n_clusters}_cluster_dtrajs.h5')
max_lag = 15

dt2 = [i.astype(np.int_) for i in dtrajs]
dt3 = [i.reshape((i.shape[0])) for i in dt2]

msm1 = msm.estimate_markov_model(dt3, 5)

print(f'Acitve state percentage is {msm1.active_state_fraction}')
print(f'Acitve count percentage is {msm1.active_count_fraction}')
