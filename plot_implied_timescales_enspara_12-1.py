import pyemma.coordinates as coor
import numpy as np
import enspara.msm as msm
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt

sys = 'fdis'
n_clusters = 100
dtrajs = coor.load(f'cluster_data/{sys}_{n_clusters}_cluster_dtrajs.h5')
max_lag = 15

dt2 = [i.astype(np.int_) for i in dtrajs]
dt3 = [i.reshape((i.shape[0])) for i in dt2]

lags = [2,4,6,8,10,12,16,24,32,48,64,80]

def norm_pseudo(C, prior_counts=1, calculate_eq_probs=True):
    return msm.builders.normalize(C, prior_counts=prior_counts, calculate_eq_probs=calculate_eq_probs)

its1 = msm.timescales.implied_timescales(np.array(dt3), lags, msm.builders.normalize, n_times=8)

fig, ax = plt.subplots()
for i in range(8):
    ax.plot(lags, np.absolute(its1[:,i]), linewidth=2, marker='o')
ax.plot(lags,lags, linestyle='dashed', linewidth=2, color='k')
ax.set_yscale('log')
ax.set_xlabel('Lag time (ns)')
ax.set_ylabel('Implied Timescales (ns)')
fig.savefig(f'{sys}_implied_timescale_enspara_norm_more3.pdf')
