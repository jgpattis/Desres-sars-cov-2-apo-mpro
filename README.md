This repository will analyze the DESRES SARS-CoV-2 protease data

scripts should be run in this order:

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
