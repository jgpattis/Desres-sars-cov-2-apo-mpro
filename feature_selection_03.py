#! /usr/bin/env/ python
# Featurize trajectories several different ways to test in next step

import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import os
import pandas as pd

feature_list = pd.read_pickle('feature_list_1.pickl')

def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):
    """Compute a cross-validated VAMP2 score.

    We randomly split the list of independent trajectories into
    a training and a validation set, compute the VAMP2 score,
    and repeat this process several times.

    Parameters
    ----------
    data : list of numpy.ndarrays
        The input data.
    dim : int
        Number of processes to score; equivalent to the dimension
        after projecting the data with VAMP2.
    lag : int
        Lag time for the VAMP2 scoring.
    number_of_splits : int, optional, default=10
        How often do we repeat the splitting and score calculation.
    validation_fraction : int, optional, default=0.5
        Fraction of trajectories which should go into the validation
        set during a split.
    """
    # we temporarily suppress very short-lived progress bars
    with pyemma.util.contexts.settings(show_progress_bars=False):
        if type(data) == list:
            nval = int(len(data) * validation_fraction)
        elif data._is_reader == True:
            nval = data.number_of_trajectories()
        else:
            raise ValueError("data must be list of numpy arrays or pyemma reader object")
        scores = np.zeros(number_of_splits)
        for n in range(number_of_splits):
            if type(data) == list:
                ival = np.random.choice(len(data), size=nval, replace=False)
            elif data._is_reader == True:
                ival = np.random.choice(data.number_of_trajectories(), size=nval, replace=False)
            vamp = coor.vamp(
                [d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
            scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
        return scores

def run_vamp_score(feat_option_list, lag_list=[2, 8, 16], dim=5, number_of_splits=10):
    '''run vamp score on list of feature options
    at several different lag times
    returns scores and errors'''
    scores = [0] * len(lag_list)
    errors = [0] * len(lag_list)
    for i, lag in enumerate(lag_list):
        scores[i] = []
        errors[i] = []
        vamp_score = [0] * len(feat_option_list)
        for j in range(len(feat_option_list)):
            vamp_score[j] = score_cv(feat_option_list[j], lag=lag, dim=dim, number_of_splits=number_of_splits)
            scores[i] += [vamp_score[j].mean()]
            errors[i] += [vamp_score[j].std()]
    return scores, errors

feat_optionsA = [coor.load(i) for i in feature_list['path_a']]
feat_optionsB = [coor.load(i) for i in feature_list['path_b']]
feat_options = [0] * len(feat_optionsA)
for i in range(len(feat_optionsA)):
    feat_options[i] = feat_optionsA[i] + feat_optionsB[i]

scores, errors = run_vamp_score(feat_options, dim=10, number_of_splits=20)

with open('vamp_scores_10dim_both_1.npy','wb') as handle:
    np.save(handle, scores)

with open('vamp_errors_10dim_both_1.npy','wb') as handle:
    np.save(handle, errors)

print(scores)
print(errors)
