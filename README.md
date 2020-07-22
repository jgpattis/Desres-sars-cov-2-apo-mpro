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
         -> histogram_00/
    filter_distance_01.py
    filter_distance_01-1.py
        |
        -> filtered_distance_featurization_01/
    featurize_02.py
        |
         -> feature_data_02/
    feature_selection_03.py
    plot_vamp_score_04.py
    tica_05.py
        |
         -> tica_data_05/
    tica_plot_06.py
        |
         -> tica_plots_06/
    mini_cluster_07.py
    pull_out_structure_08.py
        |
         -> structure_08/
    find_num_clusters_09.py
    cluster_kmeans_10.py
        |
        -> cluster_data_10/
    check_ergodic_11.py
    plot_implied_timescales_pyemma_12.py
    plot_implied_timescales_enspara_12-1.py
    plot_impolied_timescales_msmb_12-2.py
        |
        -> implied_timescales-12/

    allostery_analysis/
    cryptic_pockets/

## This analysis was performed with:
Python 3.7.4  
pyemma 2.5.7  
mdtraj 1.9.3

#### Allostery analysis uses:  
mdentropy 0.4.0.dev0  
enspara 0.1.0

#### cryptic pockets uses:
enspara 0.1.0
