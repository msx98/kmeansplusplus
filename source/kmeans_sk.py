
import numpy as np

# sklearn kmeans
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import contingency_matrix

# pyclustering kmeans
from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import distance_metric
from pyclustering.cluster.center_initializer import random_center_initializer
from pyclustering.cluster.encoder import type_encoding
from pyclustering.cluster.encoder import cluster_encoder

from typing import List


def KmeansAlgorithm(
            initial_centroids_list: List[List[float]],
            data: List[List[float]],
            dims_count: int,
            k: int,
            point_count: int,
            max_iter: int,
            eps: float
        ) -> List[List[float]]:
    initial_centroids_list = np.array([data[idx] for idx in initial_centroids_list])
    return np.around(KMeans(
            n_clusters=k,
            init=initial_centroids_list,
            tol=eps,
            n_init=1,
            #algorithm="full",
            random_state=0,
            max_iter=max_iter
        ).fit(data).cluster_centers_, 4)