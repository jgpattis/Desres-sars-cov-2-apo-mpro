#! /usr/bin/env/ python
# Plot tica data

import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import os
import pandas as pd
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt

sys = 'tica_plots_06/chi2/chi2'   ## What tica data do you want to look at
tica_data = coor.load('tica_data_05/chi2_tica_data.h5')  ## Where is that tica data?
tica_data_cat = np.concatenate(tica_data)

def kin_var_plot():
    ''' Plot kinetic variance vs time to identify how many tica
       components to use '''
    cumvar = np.load(open("tica_data_05/chi2_cumvar.npy","rb"))
    fig, ax = plt.subplots()
    index = range(1, len(cumvar) + 1)
    ax.plot(index, cumvar)
    ax.axhline(y=0.95, c='y')
    ax.axhline(y=0.90, c='C1')
    ax.axhline(y=0.85, c='r')
    ax.set_xlabel('Tica index', fontsize=16)
    ax.set_ylabel('Cumulative Variance (%)', fontsize=16)
    fig.tight_layout()
    return fig, ax

def heatmap(tica_data_cat=tica_data_cat, ic2=1):
    ''' Plot a heatmap of the density of the tica data '''
    fig, ax = plt.subplots()
    fig, ax, cb = pyemma_plots.plot_density(tica_data_cat[:,0], tica_data_cat[:,ic2], ax=ax, logscale=True)
    ax.scatter(tica_data_cat[0,0], tica_data_cat[0,ic2], marker='x', s=200, c='k')
    ax.set_xlabel(f'IC 1', fontsize=16)
    ax.set_ylabel(f'IC {ic2 + 1}', fontsize=16)
    cb['cbar'].set_label('Density', fontsize=16)
    fig.tight_layout()
    return fig, ax

def colorTIME(tica_data_cat=tica_data_cat, ic2=1, stride=2):
    ''' plots a scatter plot of tica data and colors by simulation time '''
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=[12.8, 4.8])
    time = list(range(100000))
    time_us = [i/1000 for i in time]
    mid = int(len(tica_data_cat)/2)
    print(mid)
    cbb = ax1.scatter(tica_data_cat[:mid:stride,0], tica_data_cat[:mid:stride,ic2], c=time_us[::stride])
    ax1.scatter(tica_data_cat[0,0], tica_data_cat[0,ic2], marker='x', s=200, c='k')
    ax1.set_xlabel(f'IC 1', fontsize=16)
    ax1.set_ylabel(f'IC {ic2 + 1}', fontsize=16)
    ax1.set_title('Chain A', fontsize=18)
    cb = ax2.scatter(tica_data_cat[mid::stride,0], tica_data_cat[mid::stride,ic2], c=time_us[::stride])
    ax2.scatter(tica_data_cat[0,0], tica_data_cat[0,ic2], marker='x', s=200, c='k')
    ax2.set_xlabel(f'IC 1', fontsize=16)
    ax2.set_title('Chain B', fontsize=18)
    cb2 = fig.colorbar(cbb, ax=ax1)
    cb2.set_label(f'Time (\u03BCs)', fontsize=16)
    cb1 = fig.colorbar(cb, ax=ax2)
    cb1.set_label(f'Time (\u03BCs)', fontsize=16)
    fig.tight_layout()
    return fig, (ax1, ax2)

def ic_vs_time(tica_data_cat=tica_data_cat, ic=0):
    ''' track the progression of an independent component over time '''
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
    time = list(range(100000))
    time_us = [i/1000 for i in time]
    mid = int(len(tica_data_cat)/2)
    ax1.plot(time_us, tica_data_cat[:mid,ic])
    ax1.set_title('Chain A', fontsize=18)
    ax1.set_ylabel(f'IC {ic + 1}', fontsize=16)
    ax2.plot(time_us, tica_data_cat[mid:,ic])
    ax2.set_title('Chain B', fontsize=18)
    ax2.set_ylabel(f'IC {ic + 1}', fontsize=16)
    ax2.set_xlabel(f'Time (\u03BCs)', fontsize=16)
    fig.tight_layout()
    return fig, (ax1, ax2)

fig1, ax1 = kin_var_plot()
fig1.savefig(sys + 'kinetic_variance_plot.pdf')

fig2, ax2 = heatmap(tica_data_cat=tica_data_cat, ic2=1)
fig2.savefig(sys + 'heatmap_ic2.pdf')

fig3, ax3 = heatmap(tica_data_cat=tica_data_cat, ic2=2)
fig3.savefig(sys + 'heatmap_ic3.pdf')

fig4, ax4 = heatmap(tica_data_cat=tica_data_cat, ic2=3)
fig4.savefig(sys + 'heatmap_ic4.pdf')

fig5, (ax51, ax52) = colorTIME(tica_data_cat=tica_data_cat, ic2=1)
fig5.savefig(sys + 'colortime_ic2.pdf')

fig6, (ax61, ax62) = colorTIME(tica_data_cat=tica_data_cat, ic2=2)
fig6.savefig(sys + 'colortime_ic3.pdf')

fig7, (ax71, ax72) = colorTIME(tica_data_cat=tica_data_cat, ic2=3)
fig7.savefig(sys + 'colortime_ic4.pdf')

fig8, (ax81, ax82) = ic_vs_time(tica_data_cat=tica_data_cat, ic=0)
fig8.savefig(sys + 'ICvsTIME_ic1.pdf')

fig9, (ax91, ax92) = ic_vs_time(tica_data_cat=tica_data_cat, ic=1)
fig9.savefig(sys + 'ICvsTIME_ic2.pdf')

fig10, (ax101, ax102) = ic_vs_time(tica_data_cat=tica_data_cat, ic=2)
fig10.savefig(sys + 'ICvsTIME_ic3.pdf')

fig11, (ax111, ax112) = ic_vs_time(tica_data_cat=tica_data_cat, ic=3)
fig11.savefig(sys + 'ICvsTIME_ic4.pdf')

fig12, (ax121, ax122) = ic_vs_time(tica_data_cat=tica_data_cat, ic=4)
fig12.savefig(sys + 'ICvsTIME_ic5.pdf')
