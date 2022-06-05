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
    print(f"Inner sys argv: {sys.argv}")
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
        print(args)

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


def KmeansPlusPlus(K: int, data: np.array) -> List[int]:
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


def extract_fit_params(args_override=None):
    k, max_iter, eps, file_name_1, file_name_2 = args_override or get_args()
    if not (validate_input_files(file_name_1, file_name_2)):
        print(MSG_ERR_INVALID_INPUT)
        raise Exception()
    datapoints_list = _read_data_as_np(file_name_1, file_name_2)
    verify_data(datapoints_list)
    initial_centroids_list = KmeansPlusPlus(k, datapoints_list)
    datapoints_list = sorted(datapoints_list, key=lambda x: float(x[0]))
    datapoints_list = [x[1:] for x in datapoints_list]
    point_count = len(datapoints_list)
    dims_count = len(datapoints_list[0])
    return (
        initial_centroids_list,
        datapoints_list,
        dims_count,
        k,
        point_count,
        max_iter,
        eps
    )


def main():
    fit_params = extract_fit_params()
    initial_centroids_list = fit_params[0]
    results = mykmeanssp.fit(*fit_params)
    print(','.join([str(x) for x in initial_centroids_list]))
    print('\n'.join([','.join([str(y) for y in x]) for x in results]))


if __name__ == '__main__':
    print("I am main")
    main()