#! /usr/bin/env/ python
# Plot Vamp2 score for different feature choices

import numpy as np
import pickle
import os
import pandas as pd
import matplotlib.pyplot as plt

# Load files
scores = np.load(open("vamp_scores_10dim_both_1.npy","rb"))
errors = np.load(open("vamp_errors_10dim_both_1.npy","rb"))
feature_list = pd.read_pickle('feature_list_1.pickl')

y_axe_zoom = True
y_lim = (7.9,11.05)
outfile = 'vamp_score_10dim_20splits_both_3_zoom.pdf'
def plot_vamp_score(axs, scores, errors, lags=None, labels=None, ylim=None, frame2ns=1, fs=15):
    '''Plot a bar plot for a list of vamp scores for different lag times'''
    for i, ax, lag in zip(range(len(lags)), axs.flat, lags):
        color_list = [f'C{i}' for i in range(len(scores[0]))]
        ax.bar(labels, scores[i], yerr=errors[i], color=color_list)
        ax.set_title(r'lag time $\tau$={:.1f}ns'.format(lag * frame2ns))
        if ylim is not None:
            ax.set_ylim(ylim)
    axs[0].set_ylabel('VAMP2 score', fontsize=fs)

labels = feature_list['short_label']

fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)
if y_axe_zoom == True:
    plot_vamp_score(axes, scores, errors, lags=[2, 8, 16], labels=labels, ylim=y_lim)
else:
    plot_vamp_score(axes, scores, errors, lags=[2, 8, 16], labels=labels)
fig.tight_layout()
fig.savefig(outfile)
