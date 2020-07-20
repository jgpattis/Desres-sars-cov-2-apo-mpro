## This repository will analyze the DESRES SARS-CoV-2 main protease trajectory data

The trajectory is of the SARS-CoV-2 main protease dimer, which for most of the analysis here is analysied as two independent simulations of the monomer

trajectory data can be downloaded here: http://www.deshawresearch.com/resources_sarscov2.html  
and is not included in this repository (warning trajectory is 11GB total)

this trajectory data should be cited as:

    D. E. Shaw Research, "Molecular Dynamics Simulations Related to SARS-CoV-2,"
    D. E. Shaw Research Technical Data, 2020.
    http://www.deshawresearch.com/resources_sarscov2.html

## Scripts should be run in this order:

    distance_histogram_00.py
        |
         -> histogram_00
    filter_distance_01.py
    filter_distance_01-1.py
    featurize_02.py
        |
         -> feature_data_02
    feature_selection_03.py
    plot_vamp_score_04.py
    tica_05.py
        |
         -> tica_data_05
    tica_plot_06.py
        |
         -> tica_plots_06
    mini_cluster.py
    pull_out_structure_07.py

    check_ergodic_08.py
    find_num_clusters_07.py
    plot_implied_timescales.py
    cluster_kmeans.py

    allostery_analysis

## This analysis was performed with:
Python 3.7.4  
pyemma 2.5.7  
mdtraj 1.9.3

#### Allostery analysis uses:  
mdentropy 0.4.0.dev0  
enspara 0.1.0
