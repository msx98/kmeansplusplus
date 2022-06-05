
import kmeans_pp
from typing import List
import numpy as np
from sklearn.datasets import load_iris
import sklearn.cluster
import matplotlib.pyplot as plt


def main():
    points = load_iris()['data']
    inertia = []
    for k in range(1, 11):
        result = sklearn.cluster.KMeans(
            n_clusters=k,
            init="k-means++",
            random_state=0
        ).fit(points)
        centroids = result.cluster_centers_
        inertia.append(calculate_inertia(centroids, points))
    x = [i for i in range(1, 11)]
    y = inertia
    fig = plt.figure()
    plt.plot(x,y)
    fig.savefig('elbow.png', dpi=fig.dpi)


def square_dist_from_closest_centroid(centroids: List[List[float]], point: List[float]) -> int:
    """ returns distance from centroid closest to point """
    centroids_point_diff_matrix = centroids - point
    square_dist_vector = np.sum(np.square(centroids_point_diff_matrix), axis=1)
    square_dist_of_closest_cluster = np.min(square_dist_vector)
    return square_dist_of_closest_cluster

def calculate_inertia(centroids: List[List[float]], points: List[List[float]]) -> float:
    return sum(square_dist_from_closest_centroid(centroids, point) for point in points)


if __name__ == "__main__":
    main()