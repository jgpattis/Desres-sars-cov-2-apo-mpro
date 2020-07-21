#cluster data into a small amount of clusters to later pull out structures
import pyemma.coordinates as coor
import numpy as np

sys = 'back'
tica_data = coor.load('tica_data_05/back_tica_data.h5')

n_clusters = 50

cl = coor.cluster_kmeans(tica_data, k=n_clusters, max_iter=50)

cl.save(f'{sys}_{n_clusters}_mini_cluster_object.h5', overwrite=True)

cl.write_to_hdf5(f'{sys}_{n_clusters}_cluster_dtrajs.h5')
