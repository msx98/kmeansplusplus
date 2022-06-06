import numpy as np
from typing import List, Tuple


def locate_closest_centroid(centroids: List[List[float]], point: List[float]) -> int:
    """ returns index of centroid closest to point """
    diff_matrix = centroids - point
    dist_vector = np.sum(diff_matrix**2, axis=1)
    idx_of_closest_cluster = np.argmin(dist_vector)
    return idx_of_closest_cluster

def get_closest_clusters_map(centroids: List[List[float]], points: List[List[float]]) -> List[int]:
    """ returns array. arr[j] = index of centroid closest to point j"""
    indices = []
    for i in range(len(points)):
        idx_of_closest_centroid = locate_closest_centroid(centroids, points[i])
        indices.append(idx_of_closest_centroid)
    return indices

def get_closest_clusters_map_np(centroids: List[List[float]], points: List[List[float]]) -> List[int]:
    centroids, points = np.array(centroids), np.array(points)
    k, dims = centroids.shape
    n, _ = points.shape
    agg_centroids = np.repeat(centroids[:, :, np.newaxis], n, axis=2)
    agg_points = points[:,:,np.newaxis].transpose(2,1,0)
    #agg_points=np.repeat(agg_points, k, axis=0)
    return np.argmin(np.sum(np.square(agg_centroids-agg_points),axis=1), axis=0)

def get_cluster_sets_np(centroids: List[List[float]], points: List[List[float]]) -> List[List[float]]:
    centroids, points = np.array(centroids), np.array(points)
    k, dims = centroids.shape
    n, _ = points.shape
    cluster_for_point = get_closest_clusters_map_np(centroids, points)
    S = [[] for _ in range(k)]
    for i in range(n):
        S[cluster_for_point[i]].append(points[i])
    new_centroids = [sum(x)/len(x) for x in S]
    return new_centroids



def get_closest_clusters_map(centroids: List[List[float]], points: List[List[float]]) -> List[int]:
    """ returns array. arr[j] = index of centroid closest to point j"""
    k = len(centroids)
    repeated_centroids = np.repeat(np.array(centroids)[:, :, np.newaxis], k, axis=2)
    indices = []
    for i in range(len(points)):
        idx_of_closest_centroid = locate_closest_centroid(centroids, points[i])
        indices.append(idx_of_closest_centroid)
    return indices

def get_cluster_sets(centroids: List[List[float]], points: List[List[float]]) -> List[List[int]]:
    """ returns [S_1, ..., S_k] """
    k = len(centroids)
    clusters = [[] for _ in range(k)]
    closest_clusters_map = get_closest_clusters_map(centroids, points)
    for point, cluster in enumerate(closest_clusters_map):
        clusters[cluster].append(point)
    for cluster in clusters:
        assert(len(cluster) > 0)
    return clusters


def update_centroids(centroids: List[List[float]], points: List[List[float]]) -> List[List[float]]:
    """ reassigns centroids """
    k = len(centroids)
    dims = len(centroids[0])
    new_centroids = []
    points_for_cluster = get_cluster_sets(centroids, points)
    S = [[points[i] for i in points_list] for points_list in points_for_cluster]
    for j in range(k):
        new_centroid = np.zeros(dims)
        for point_idx in S[j]:
            new_centroid += np.array(points[point_idx])
        new_centroid /= len(S[j])
        new_centroids.append(new_centroid)
    return new_centroids

update_centroids = get_cluster_sets_np

def is_convergence(previous_centroids: List[List[float]], centroids: List[List[float]], epsilon: float) -> bool:
    previous_centroids = np.array(previous_centroids)
    centroids = np.array(centroids)
    diff_matrix = previous_centroids-centroids
    dist_vector = np.sqrt(np.sum(diff_matrix**2, axis=1))
    return np.all(dist_vector < epsilon)

def KmeansAlgorithm(
            initial_centroids_list: List[List[float]],
            data: List[List[float]],
            dims_count: int,
            k: int,
            point_count: int,
            max_iter: int,
            eps: float
        ) -> List[List[float]]:
    centroids_list = [data[x] for x in initial_centroids_list]
    for i in range(max_iter):
        previous_centroids_list = centroids_list
        centroids_list = update_centroids(centroids_list, data)
        if is_convergence(previous_centroids_list, centroids_list, eps):
            return np.around(centroids_list,4)
    return np.around(centroids_list,4)