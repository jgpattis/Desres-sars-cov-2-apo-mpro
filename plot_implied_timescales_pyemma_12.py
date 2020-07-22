import pyemma.coordinates as coor
import numpy as np
import pyemma.msm as msm
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt

sys = 'fdis'
n_clusters = 100
dtrajs = coor.load(f'cluster_data/{sys}_{n_clusters}_cluster_dtrajs.h5')
max_lag = 80

dt2 = [i.astype(np.int_) for i in dtrajs]
dt3 = [i.reshape((i.shape[0])) for i in dt2]

its = msm.its(dt3, lags=max_lag, nits=8, errors='bayes', nsamples=200)

fig, ax = plt.subplots()
pyemma_plots.plot_implied_timescales(its, units='ns', ax=ax)
fig.savefig(f'{sys}_implied_timescale_{max_lag}.pdf')
