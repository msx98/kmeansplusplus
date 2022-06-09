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


def main():
    np.random.seed(0)
    fit_params = extract_fit_params()
    initial_centroids_list = fit_params[0]
    results = mykmeanssp.fit(*fit_params)
    print(','.join([str(x) for x in initial_centroids_list]))
    print('\n'.join([','.join(["%.4f"%y for y in x]) for x in results]))


def extract_fit_params():
    k, max_iter, eps, file_name_1, file_name_2 = get_cmd_args()
    if not (validate_input_files(file_name_1, file_name_2)):
        print(MSG_ERR_INVALID_INPUT)
        raise Exception()
    datapoints_list = _read_data_as_np(file_name_1, file_name_2)
    verify_data(datapoints_list)
    if k >= len(datapoints_list):
        print(MSG_ERR_INVALID_INPUT)
        raise Exception()
    initial_centroids_list = KMeansPlusPlus(k, datapoints_list)
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

def _find_first_centroids(K: int, data: np.array) -> List[int]:
    np.hstack((np.array(range(1, data.shape[0]+1)).reshape((data.shape[0],1)), data))
    number_of_points_in_data_N = len(data)  # the number of rows(points) in data
    # todo check for error from last time
    if K >= number_of_points_in_data_N:
        print(MSG_ERR_INVALID_INPUT)
        raise Exception()
    options = list(data[:, 0])
    np.random.seed(0)
    centroid_list = [np.random.choice(options)]
    
    min_distance_vector = np.zeros(number_of_points_in_data_N)
    probability_vector = np.zeros(number_of_points_in_data_N)

    for i in range(K):
        pass

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


def KMeansPlusPlus_2(k: int, x: np.array) -> List[int]:
    np.random.rand(0)
    N, d = x.shape
    u = [None for _ in range(N)]
    u_idx = [-1 for _ in range(N)]
    P = [0 for _ in range(N)]
    D = [float('inf') for _ in range(N)]

    i = 0
    selection = int(np.random.choice(x[:,0]))
    u[0] = x[selection]

    while i < k:
        for l in range(N):
            x_l = x[l] # remove index
            min_square_dist = float('inf')
            for j in range(0,i+1):
                u_j = u[j] # remove index
                square_dist = np.sum((x_l[1:] - u_j[1:])**2)
                min_square_dist = min(square_dist, min_square_dist)
            D[l] = min_square_dist
        D_sum = sum(D)
        P = D/D_sum

        i += 1
        selection = int(np.random.choice(x[:,0], p=P))
        u[i] = x[selection]
        continue

    return u


def KMeansPlusPlus(k: int, data: np.array) -> List[int]:
    return KMeansPlusPlus_2(k, data)
    #return _find_first_centroids(k, data)
    data = np.copy(data)
    N, dims = data.shape

    D = np.ones(N)*np.inf
    P = np.ones(N)/N
    centroids_list = np.ones((k, dims)) * np.inf
    centroids_choice = np.ones(k, dtype=np.uint64) * (-1)

    np.random.seed(0)

    for i in range(1,k):
        centroids_choice[i] = np.random.choice(N, p=P)
        centroids_list[i] = data[int(centroids_choice[i])]
        for l in range(N):
            D[l] = np.min(np.sum(np.square(data[l]-centroids_list), axis=1))
        P = D/np.sum(D)
    return centroids_choice


def _read_data_as_np(file_name1: str, file_name2: str) -> np.array:
    path_file1 = os.path.join(os.getcwd(), file_name1)
    path_file2 = os.path.join(os.getcwd(), file_name2)
    data_frame_1 = pd.read_csv(path_file1, header=None).rename({0: "index"}, axis=1)
    data_frame_2 = pd.read_csv(path_file2, header=None).rename({0: "index"}, axis=1)
    joined_data_frame = data_frame_1.join(
                                        data_frame_2.set_index('index'),
                                        on='index', lsuffix='from_second_file ',
                                        how='inner')
    joined_data_frame = joined_data_frame.sort_values('index')
    #joined_data_frame.drop('index', inplace=True, axis=1)
    data = joined_data_frame.to_numpy()
    return data


def validate_input_files(file_name1: str, file_name2: str) -> bool:
    for file in [file_name1, file_name2]:
        if not (file.endswith("csv") or file.endswith("txt")):
            return False
    return True


def verify_data(data: List[List[float]]):
    N = len(data)
    if len(data) == 0:
        print(MSG_ERR_GENERIC)
        exit(1)
    dims = len(data[0])
    for point in data:
        if len(point) != dims:
            print(MSG_ERR_GENERIC)
            exit(1)


def get_cmd_args():
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
        print(args)

        print(MSG_ERR_INVALID_INPUT)
        exit(1)


if __name__ == '__main__':
    main()
