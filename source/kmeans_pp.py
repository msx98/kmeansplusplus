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
    fit_params = extract_fit_params()
    initial_centroids_list = fit_params[0]
    results = mykmeanssp.fit(*fit_params)
    print(','.join([str(x) for x in initial_centroids_list]))
    print('\n'.join([','.join([str(y) for y in x]) for x in results]))


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


def KMeansPlusPlus(k: int, data: np.array) -> List[int]:
    data = np.copy(data)
    N = len(data)
    dims = len(data[0])

    D = np.ones(N)*np.inf
    P = np.ones(N)/N
    centroids_list = np.ones((k, dims))*np.inf
    centroids_choice = np.ones(k, dtype=np.uint64)

    for i in range(k):
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
                                        how='left')
    joined_data_frame = joined_data_frame.sort_values('index')
    joined_data_frame.drop('index', inplace=True, axis=1)
    data = joined_data_frame.to_numpy()
    return data


def validate_input_files(file_name1: str, file_name2: str) -> bool:
    for file in [file_name1, file_name2]:
        if not (file.endswith("csv") or file.endswith("txt")):
            return False
    return True


def verify_data(data: List[List[float]]):
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