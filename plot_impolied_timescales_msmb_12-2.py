import pyemma.coordinates as coor
import numpy as np
from multiprocessing import Pool
import pandas as pd
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt
from msmbuilder.msm import MarkovStateModel
import seaborn as sns

sys = 'fdis'
n_clusters = 100
dtrajs = coor.load(f'cluster_data/{sys}_{n_clusters}_cluster_dtrajs.h5')
max_lag = 80

dt2 = [i.astype(np.int_) for i in dtrajs]
dt3 = [i.reshape((i.shape[0])) for i in dt2]

#lagtimes = [2 ** i for i in range(8)]
lagtimes = [2,4,6,8,10,12,16,24,32,48,64,80] 

## Define what to do for parallel execution
def at_lagtime(lt):
    msm = MarkovStateModel(lag_time=lt, n_timescales=8, ergodic_cutoff='off', reversible_type='none', verbose=False)
    msm.fit(dt3)
    ret = {
        'lag_time': lt,
        'percent_retained': msm.percent_retained_,
    }
    for i in range(msm.n_timescales):
        ret['timescale_{}'.format(i)] = msm.timescales_[i]
    return ret


## Do the calculation
with Pool() as p:
    results = p.map(at_lagtime, lagtimes)

timescales = pd.DataFrame(results)

n_timescales = 8
sns.set_style('ticks')
colors = sns.color_palette()

## Implied timescales vs lagtime
def plot_timescales(ax):
    for i in range(n_timescales):
        ax.plot(timescales['lag_time'],
                   np.absolute(timescales['timescale_{}'.format(i)]),
                   linewidth=2, marker='o'
                   )
    
    xmin, xmax = ax.get_xlim()
    xx = np.linspace(xmin, xmax)
    ax.plot(xx, xx, color=colors[2], label='$y=x$')
    ax.legend(loc='best', fontsize=14)
    ax.set_xlabel('Lag Time (ns)', fontsize=18)
    ax.set_ylabel('Implied Timescales (ns)', fontsize=18)
    #ax.set_xscale('log')
    ax.set_yscale('log')

## Percent trimmed vs lagtime
def plot_trimmed(ax):
    ax.plot(timescales['lag_time'],
            timescales['percent_retained'],
            'o-',
            label=None,  # pandas be interfering
            )
    ax.axhline(100, color='k', ls='--', label='100%')
    ax.legend(loc='best', fontsize=14)
    ax.set_xlabel('Lag Time (ns)', fontsize=18)
    ax.set_ylabel('Percent Retained', fontsize=18)
    #ax.set_xscale('log')
    ax.set_ylim((0, 110))

## Plot timescales
fig, ax = plt.subplots(figsize=(7, 5))
plot_timescales(ax)
fig.tight_layout()
fig.savefig('implied-timescales_msmb_off.pdf')

## Plot trimmed
fig, ax = plt.subplots(figsize=(7,5))
plot_trimmed(ax)
fig.tight_layout()
fig.savefig('percent-trimmed_msmb_off.pdf')
