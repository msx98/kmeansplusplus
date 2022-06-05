import os
import sys
from math import sqrt
from typing import List, Dict, Set
import numpy as np

MSG_ERR_INVALID_INPUT = "Invalid Input!"
MSG_ERR_GENERIC       = "An Error Has Occurred"

EPSILON = 0.001
INFINITY = float('inf')
MAX_ITER_UNSPEC = 200


def main():
    k, max_iter, path_to_input, path_to_output = get_args()
    KmeansAlgorithm(k, max_iter, path_to_input, path_to_output)


def get_args():
    args = sys.argv
    if args[0] in ["python", "python3", "python.exe", "python3.exe"]:
        args = args[1:]
    if args[0][-3:] == ".py":
        args = args[1:]
    try:
        if len(args) == 3:
            return int(args[0]), MAX_ITER_UNSPEC, args[1], args[2]
        elif len(args) == 4:
            return int(args[0]), int(args[1]), args[2], args[3]
        else:
            raise Exception()
    except:
        print(MSG_ERR_INVALID_INPUT)
        exit(1)


def KmeansAlgorithm_Files(K: int, max_iter=200, input_filename: str = None, output_filename: str = None) -> str:
    file_lines = None
    data = _read_data(input_filename)
    verify_data(data)

    centroids_list = data[0:K]

    for i in range(max_iter):
        point_to_centroid_list = _point_to_centroid_list(data, centroids_list)
        previous_centroids_list = centroids_list
        centroids_list = _update_centroid(centroids_list, data, point_to_centroid_list)
        if _is_convergence(previous_centroids_list, centroids_list):
            # print("reach convergence")
            return _write_centroid_to_text(output_filename, centroids_list)
    return _write_centroid_to_text(output_filename, centroids_list)


def KmeansAlgorithm(
            initial_centroids_list: List[List[float]],
            data: List[List[float]],
            dims_count: int,
            k: int,
            point_count: int,
            max_iter: int,
            eps: float
        ) -> List[List[float]]:
    verify_data(data)
    centroids_list = [data[x] for x in initial_centroids_list]
    for i in range(max_iter):
        point_to_centroid_list = _point_to_centroid_list(data, centroids_list)
        previous_centroids_list = centroids_list
        centroids_list = _update_centroid(centroids_list, data, point_to_centroid_list)
        if _is_convergence(previous_centroids_list, centroids_list, eps):
            # print("reach convergence")
            return np.around(centroids_list,4)
    return np.around(centroids_list,4)


def _read_data(input_filename) -> List[List[float]]:
    # creating the path of the file that we read.
    # getcwd return the path of the current directory
    path = os.path.join(os.getcwd(), input_filename)
    try:
        with open(path, 'r') as file:
            file_lines = file.readlines()
            # vector_lst -> [vector[string,...,string],...,]]
            vector_lst = [temp_line.split(',') for temp_line in file_lines]
            # data -> list[vector[float,...,float],..,]
            data = [list(map(float, vector)) for vector in vector_lst]
        return data
    except FileNotFoundError as e:
        print(MSG_ERR_GENERIC)
        raise FileNotFoundError()

# this function verifies that all points have the same dimension
# and also that we have a nonzero number of points
def verify_data(data: List[List[float]]):
    if len(data) == 0:
        print(MSG_ERR_GENERIC)
        exit(1)
    dims = len(data[0])
    for point in data:
        if len(point) != dims:
            print(MSG_ERR_GENERIC)
            exit(1)


# A function that calculate the distance between 2 points.

def _distance_between_point_and_centroid(point: List[float], centroid: List[float]) -> float:
    distance = 0
    for i in range(len(point)):
        distance += (point[i] - centroid[i]) ** 2
    return distance ** 0.5


def _point_to_centroid_list(data: List[List[float]], centroids_list: List[List[float]]) -> List[List[float]]:
    mapping_centroid_to_each_point = []
    for index, point in enumerate(data):
        minimum_distance = sys.float_info.max  # max value in python
        curr_centroid = None
        for centroid in centroids_list:
            curr_distance = _distance_between_point_and_centroid(point, centroid)
            if curr_distance < minimum_distance:
                minimum_distance = curr_distance
                curr_centroid = centroid
        mapping_centroid_to_each_point.append(curr_centroid)

    return mapping_centroid_to_each_point


def _update_centroid(centroids: List[List[float]], data: List[List[float]],
                     mapping_point_to_centroid: List[List[float]]) -> List[List[float]]:
    updated_centroid_list = []
    for centroid in centroids:
        new_updated_centroid = []
        # creating the centroid cluster
        centroid_cluster = [data[i] for i in range(len(data)) if np.all(mapping_point_to_centroid[i] == centroid)]
        cluster_size = len(centroid_cluster)
        for i in zip(*centroid_cluster):
            cluster_sum = sum(i)
            new_updated_centroid.append(cluster_sum / cluster_size)
        updated_centroid_list.append(new_updated_centroid)

    return updated_centroid_list


def _is_convergence(prev_centroids: List[List[float]], curr_centroids: List[List[float]], eps:float=EPSILON) -> bool:
    for i in range(len(prev_centroids)):
        dis_between_2_centroids = _distance_between_point_and_centroid(prev_centroids[i], curr_centroids[i])
        if dis_between_2_centroids > eps:
            return False
    # if the difference between all the centroids have changed less than epsilon -> stop.
    return True


def _write_centroid_to_text(output_file: str, centroid_list: List[List[float]]):
    output_path = os.path.join(os.getcwd(), output_file)
    with open(output_path, 'w') as file:
        file.writelines(','.join(_format_point(point) for point in centroid) + '\n' for centroid in centroid_list)


def _format_point(point):
    return "%.4f" % point


#try: os.mkdir("resources/out")
#except: pass
#KmeansAlgorithm_Files(3, 600, "resources/input_1.txt", "resources/out/output_1.txt")
#KmeansAlgorithm_Files(7, 200, "resources/input_2.txt", "resources/out/output_2.txt")
#KmeansAlgorithm_Files(15, 300, "resources/input_3.txt", "resources/out/output_3.txt")

if __name__ == "main":
    main()