
import kmeans_pp
from typing import List
import numpy as np
from sklearn.datasets import load_iris
import sklearn.cluster
import matplotlib.pyplot as plt
import math


def main():
    points = load_iris()['data']
    inertia = []
    for k in range(2, 11):
        result = sklearn.cluster.KMeans(
            n_clusters=k,
            init="k-means++",
            random_state=0
        ).fit(points)
        centroids = result.cluster_centers_
        inertia.append(calculate_inertia(centroids, points))
    k_elbow = find_elbow(inertia) + 2
    k_elbow_inertia = inertia[k_elbow-1]
    print(k_elbow)
    x = [i for i in range(2, 11)]
    y = inertia
    fig = plt.figure()
    plt.plot(x,y)
    plt.scatter(k_elbow, k_elbow_inertia, s=80, facecolors='none', edgecolors='black', linestyle='--')
    plt.annotate("Elbow Point", xy=(k_elbow, k_elbow_inertia), xytext=(k_elbow+1, k_elbow_inertia),
            arrowprops=dict(arrowstyle="->", linestyle="--"))
    plt.xlabel("k")
    plt.ylabel("Average Dispersion")
    plt.title("Elbow Method for selection of optimal \"K\" clusters")
    fig.savefig('elbow.png', dpi=fig.dpi)


def find_elbow(x: List[float]) -> int:
    max_second_derivative_idx = 0
    max_second_derivative = -1 * math.inf
    for i in range(1, len(x)-1):
        second_derivative = calculate_second_derivative(x, i)
        if second_derivative > max_second_derivative:
            max_second_derivative_idx = i
            max_second_derivative = second_derivative
    return max_second_derivative_idx

def square_dist_from_closest_centroid(centroids: List[List[float]], point: List[float]) -> int:
    """ returns distance from centroid closest to point """
    centroids_point_diff_matrix = centroids - point
    square_dist_vector = np.sum(np.square(centroids_point_diff_matrix), axis=1)
    square_dist_of_closest_cluster = np.min(square_dist_vector)
    return square_dist_of_closest_cluster

def calculate_inertia(centroids: List[List[float]], points: List[List[float]]) -> float:
    return sum(square_dist_from_closest_centroid(centroids, point) for point in points)

def calculate_second_derivative(x: List[float], i: int) -> float:
    return x[i+1] + x[i-1] - 2 * x[i]


if __name__ == "__main__":
    main()
