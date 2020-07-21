#! /usr/bin/env/ python
# Featurize trajectories several different ways to test in next step

import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import os
import pandas as pd
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt

sys = 'fdis'
tica_data = coor.load('tica_data_05/fdis_tica_data.h5')

n_clustercenters = [25, 50, 100, 200, 400, 600]#, 800, 1000, 1200]

scores = np.zeros((len(n_clustercenters), 10))
for n, k in enumerate(n_clustercenters):
    for m in range(10):
        with pyemma.util.contexts.settings(show_progress_bars=False):
            _cl = pyemma.coordinates.cluster_mini_batch_kmeans(
                tica_data, k=k, max_iter=50)
            _msm = pyemma.msm.estimate_markov_model(_cl.dtrajs, 5)
            scores[n, m] = _msm.score_cv(
                _cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k))

with open('scores_cluster_num2.pkl','wb') as handle:
        pickle.dump(scores, handle)

fig, ax = plt.subplots()
lower, upper = pyemma.util.statistics.confidence_interval(scores.T.tolist(), conf=0.9)
ax.fill_between(n_clustercenters, lower, upper, alpha=0.3)
ax.plot(n_clustercenters, np.mean(scores, axis=1), '-o')
#ax.semilogx()
ax.set_xlabel('number of cluster centers')
ax.set_ylabel('VAMP-2 score')
fig.tight_layout()
fig.savefig(f'{sys}_cluster_vamp_score-3.pdf')
