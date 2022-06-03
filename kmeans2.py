import os
import sys
from typing import List

import numpy as np
import pandas as pd

MSG_ERR_INVALID_INPUT = "Invalid Input!"
MSG_ERR_GENERIC = "An Error Has Occurred"

EPSILON = 0.001
INFINITY = float('inf')
MAX_ITER_UNSPEC = 300


def get_args():
    args = sys.argv[2:]
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


def K_Mean_Plus():
    K, max_iter, epsilon, file_name1, file_name2 = get_args()
    if not (validate_input_files(file_name1, file_name2)):
        print(MSG_ERR_INVALID_INPUT)
        raise Exception()
    data = _read_data_as_np(file_name1, file_name2)
    first_centroids = _find_first_centroids(K, data)


def _find_first_centroids(K: int, data: np.array) -> List[int]:
    number_of_points_in_data_N = len(data)  # the number of rows(points) in data
    # todo check for error from last time
    if K >= number_of_points_in_data_N:
        print(MSG_ERR_INVALID_INPUT)
        raise Exception()
    np.random.seed(0)
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


if __name__ == '__main__':
    # print(_find_first_centroids(2))
    data = _read_data_as_np("input_1_db_1.txt", "input_1_db_2.txt")
    print(_find_first_centroids(3, data))

