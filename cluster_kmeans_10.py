#cluster tica data into clusters
import pyemma.coordinates as coor
import numpy as np

sys = 'fdis'
tica_data = coor.load('tica_data_05/fdis_tica_data.h5')

n_clusters = 100

cl = coor.cluster_kmeans(tica_data, k=n_clusters, max_iter=50)

#cl.save(f'cluster_data/{sys}_{n_clusters}_mini_cluster_object.h5', overwrite=True)

cl.write_to_hdf5(f'cluster_data_11/{sys}_{n_clusters}_cluster_dtrajs22.h5')
