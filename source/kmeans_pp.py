import os
import sys
from typing import List
import numpy as np
import pandas as pd
import mykmeanssp

MSG_ERR_INVALID_INPUT = "Invalid Input!"
MSG_ERR_GENERIC = "An Error Has Occurred"

EPSILON = 0.001
INFINITY = float('inf')
MAX_ITER_UNSPEC = 300

np.random.seed(0)

"""
k = 10
dims_count = 3
point_count = 25
max_iter=10
eps=0.01

datapoints_list = np.random.rand(point_count, dims_count)
initial_centroids_list = np.random.rand(k, dims_count)

result = mykmeanssp.fit(
    initial_centroids_list,
    datapoints_list,
    dims_count,
    k,
    point_count,
    max_iter,
    eps
)
"""

def get_args():
    args = sys.argv
    if args[0] in ["python", "python3", "python.exe", "python3.exe"]:
        args = args[1:]
    if args[0][-3:] == ".py":
        args = args[1:]
    try:
        if len(args) == 4:  # without max_itr
            return int(args[0]), MAX_ITER_UNSPEC, float(args[1]), args[2], args[3]
        elif len(args) == 5:
            return int(args[0]), int(args[1]), float(args[2]), args[3], args[4]
        else:
            raise Exception()
    except:
        print(MSG_ERR_INVALID_INPUT)
        exit(1)


def validate_input_files(file_name1: str, file_name2: str) -> bool:
    for file in [file_name1, file_name2]:
        if not (file.endswith("csv") or file.endswith("txt")):
            return False
    return True


def _read_data_as_np(file_name1: str, file_name2: str) -> np.array:
    # creating the path of the file that we read.
    # getcwd return the path of the current directory
    path_file1 = os.path.join(os.getcwd(), file_name1)
    path_file2 = os.path.join(os.getcwd(), file_name2)

    data_frame_1 = pd.read_csv(path_file1, header=None).rename({0: "index"}, axis=1)
    data_frame_2 = pd.read_csv(path_file2, header=None).rename({0: "index"}, axis=1)

    joined_data_frame = data_frame_1.join(data_frame_2.set_index('index'), on='index', lsuffix='from_second_file ',
                                          how='left')
    data = joined_data_frame.to_numpy()
    return data


def _find_first_centroids(K: int, data: np.array) -> List[int]:
    number_of_points_in_data_N = len(data)  # the number of rows(points) in data
    # todo check for error from last time
    if K >= number_of_points_in_data_N:
        print(MSG_ERR_INVALID_INPUT)
        raise Exception()
    centroid_list = [int(data[0][0])]

    min_distance_vector = np.zeros(number_of_points_in_data_N)
    probability_vector = np.zeros(number_of_points_in_data_N)

    for i in range(1, K):
        for point in data:
            minimum_distance = np.inf
            for index in centroid_list:
                centroid_index = np.where(data[:, 0] == index)
                centroid = data[centroid_index][0]
                current_distance = np.sqrt(np.sum((point[1:] - centroid[1:]) ** 2))
                minimum_distance = min(current_distance, minimum_distance)
            min_distance_vector[int(point[0])] = minimum_distance
        mini_distance_sum = sum(min_distance_vector)  # sum all args in min_disance_vector (denominator)
        probability_vector = min_distance_vector / mini_distance_sum  # calculate the probability for each
        # todo np.max or this calcualte (choose random or the maximum) 
        random_prob_choice = np.random.choice(probability_vector)
        next_centroid_index = int(np.where(probability_vector == random_prob_choice)[0])
        centroid_list.append(next_centroid_index)

    return centroid_list

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

# A function that calculates the distance between 2 points.
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
            new_updated_centroid.append(round(cluster_sum / cluster_size, 4))
        updated_centroid_list.append(new_updated_centroid)

    return updated_centroid_list


def _is_convergence(prev_centroids: List[List[float]], curr_centroids: List[List[float]], eps: float) -> bool:
    for i in range(len(prev_centroids)):
        dis_between_2_centroids = _distance_between_point_and_centroid(prev_centroids[i], curr_centroids[i])
        if dis_between_2_centroids > eps:
            return False
    # if the difference between all the centroids have changed less than epsilon -> stop.
    return True

def KmeanAlgorithm(
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
            return centroids_list
    return centroids_list

def KmeanAlgorithm_C(
            initial_centroids_list: List[List[float]],
            data: List[List[float]],
            dims_count: int,
            k: int,
            point_count: int,
            max_iter: int,
            eps: float
        ) -> List[List[float]]:
    return mykmeanssp.fit(
        initial_centroids_list,
        datapoints_list,
        dims_count,
        k,
        point_count,
        max_iter,
        eps
    )


if __name__ == '__main__':
    # print(_find_first_centroids(2))
    k, max_iter, eps, file_name_1, file_name_2 = get_args()
    print(f"Set args: k={k}, max_iter={max_iter}, eps={eps}")
    datapoints_list = _read_data_as_np(file_name_1, file_name_2)
    verify_data(datapoints_list)
    initial_centroids_list = _find_first_centroids(k, datapoints_list)
    print(f"Initial: {initial_centroids_list}")


    desired = None
    print("Expected:")
    file_expected = file_name_1.split("_")[0].replace("input","output")+"_"+file_name_1.split("_")[1]+".txt"
    with open(file_expected, 'r') as f:
        s = f.read().split("\n")[:-1]
        s = [x.split(",") for x in s]
        s[0] = [int(y) for y in s[0]]
        initial_centroids_list = s[0]
        print("Desired: " + str([int(y) for y in s[0]]))
        desired = [[float(y) for y in x] for x in s[1:]]
    print(desired)
    
    print(f"set initial centroids: {initial_centroids_list}")

    datapoints_list = sorted(datapoints_list, key=lambda x: float(x[0]))
    datapoints_list = [x[1:] for x in datapoints_list]

    point_count = len(datapoints_list)
    dims_count = len(datapoints_list[0])

    centroids_list = KmeanAlgorithm(
        initial_centroids_list,
        datapoints_list,
        dims_count,
        k,
        point_count,
        max_iter,
        eps
    )
    print("Py output:")
    print(centroids_list)

    centroids_list_c = KmeanAlgorithm_C(
        initial_centroids_list,
        datapoints_list,
        dims_count,
        k,
        point_count,
        max_iter,
        eps
    )
    centroids_list_c = [[y for y in x] for x in centroids_list_c]
    print("C output:")
    print(centroids_list_c)

    desired = np.array(desired)
    centroids_list = np.array(centroids_list)
    centroids_list_c = np.array(centroids_list_c)

    dist_py = np.all(np.abs(desired-centroids_list) < 0.001)
    dist_c  = np.all(np.abs(desired-centroids_list_c) < 0.001)
    assert(dist_py and dist_c)
