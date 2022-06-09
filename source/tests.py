#!/usr/bin/env python3

import unittest
from unittest.mock import patch
import pickle
import logging
import sys
import numpy as np
import sklearn.cluster
import kmeans_pp
import kmeans as kmeans_2
import kmeans_mine as kmeans_1
import kmeans_sk
import mykmeanssp
from typing import List
import math
import time
from sklearn.datasets import make_blobs
import re


class TestFit(unittest.TestCase):

    @unittest.skip("blah")
    def test_c_and_sklearn_equal(self):
        print("test_py_and_c_equal_random() - start")
        seed = np.random.randint(1, 100000)
        print(f"setting seed={seed}")
        np.random.seed(seed)
        fit_params = list(randomize_fit_params(k=10, max_iter=150, eps=0.01))
        (
            initial_centroids_list,
            datapoints_list,
            dims_count,
            k,
            point_count,
            max_iter,
            eps
        ) = fit_params
        centroids_list_py = kmeans_1.KmeansAlgorithm(*fit_params)
        centroids_list_c = mykmeanssp.fit(*fit_params)
        centroids_list_sk = sklearn.cluster.KMeans(
            n_clusters=k,
            init=np.array([datapoints_list[idx] for idx in initial_centroids_list]),
            max_iter=max_iter,
            tol=eps,
            n_init=1,
            random_state=0,
            algorithm="full",
        ).fit(datapoints_list).cluster_centers_
        centroids_list_sk = np.sort(np.array(centroids_list_sk))
        centroids_list_c = np.sort(np.array(centroids_list_c))
        centroids_list_py = np.sort(np.array(centroids_list_py))
        #print(f"centroids_list_sk={centroids_list_sk}")
        #print(f"ceeeeeee: {centroids_list_c}")
        relative_error = round((np.linalg.norm(centroids_list_sk-centroids_list_c)/np.linalg.norm(centroids_list_sk))*100,2)
        relative_error_py_c = round((np.linalg.norm(centroids_list_py-centroids_list_c)/np.linalg.norm(centroids_list_py))*100,2)
        relative_error_py_sk = round((np.linalg.norm(centroids_list_py-centroids_list_sk)/np.linalg.norm(centroids_list_sk))*100,2)
        print(f"relative errors:  c,sk = {relative_error}%,  c,py = {relative_error_py_c}%,  py,sk = {relative_error_py_sk}%")
        #print(centroids_list_sk)
        #print(centroids_list_py)
        """
        if relative_error_py_c > 0:
            print("blah")
            print(centroids_list_c)
            print(centroids_list_py)
            print(f" \n\
            (\n\
            initial_centroids_list={initial_centroids_list},\n\
            datapoints_list={datapoints_list},\n\
            dims_count={dims_count},\n\
            k={k},\n\
            point_count={point_count},\n\
            max_iter={max_iter},\n\
            eps={eps}\n\
            )")
            self.assertTrue(False)
        else:
            #print("Durkha")
            pass
        self.assertTrue(True)#relative_error < 5.00)
        """
    
    @unittest.skip("blah")
    def test_py_and_mine_equal(self):
        print("test_py_and_mine_equal_random() - start")
        seed = np.random.randint(1, 100000)
        print(f"setting seed={seed}")
        np.random.seed(seed)
        fit_params = list(randomize_fit_params(k=3, max_iter=150, eps=0.01))
        (
            initial_centroids_list,
            datapoints_list,
            dims_count,
            k,
            point_count,
            max_iter,
            eps
        ) = fit_params
        print(initial_centroids_list)
        centroids_list_py = kmeans_1.KmeansAlgorithm(*fit_params)
        centroids_list_c = kmeans_2.KmeansAlgorithm(*fit_params)
        centroids_list_sk = sklearn.cluster.KMeans(
            n_clusters=k,
            init=np.array([datapoints_list[idx] for idx in initial_centroids_list]),
            max_iter=max_iter,
            tol=eps,
            n_init=1,
            random_state=0,
            algorithm="full",
        ).fit(datapoints_list).cluster_centers_
        centroids_list_sk = np.array(centroids_list_sk)
        centroids_list_c = np.array(centroids_list_c)
        centroids_list_py = np.array(centroids_list_py)
        #print(f"centroids_list_sk={centroids_list_sk}")
        #print(f"ceeeeeee: {centroids_list_c}")
        relative_error = round((np.linalg.norm(centroids_list_sk-centroids_list_c)/np.linalg.norm(centroids_list_sk))*100,2)
        relative_error_py_c = round((np.linalg.norm(centroids_list_py-centroids_list_c)/np.linalg.norm(centroids_list_py))*100,2)
        relative_error_py_sk = round((np.linalg.norm(centroids_list_py-centroids_list_sk)/np.linalg.norm(centroids_list_sk))*100,2)
        print(f"relative errors:  mine,sk = {relative_error}%,  mine,py = {relative_error_py_c}%,  py,sk = {relative_error_py_sk}%")
        #print(centroids_list_sk)
        #print(centroids_list_py)
        """
        if relative_error_py_c > 0:
            print("blah")
            print(centroids_list_c)
            print(centroids_list_py)
            print(f" \n\
            (\n\
            initial_centroids_list={initial_centroids_list},\n\
            datapoints_list={datapoints_list},\n\
            dims_count={dims_count},\n\
            k={k},\n\
            point_count={point_count},\n\
            max_iter={max_iter},\n\
            eps={eps}\n\
            )")
            self.assertTrue(False)
        else:
            #print("Durkha")
            pass"""
        self.assertTrue(relative_error_py_c == 0)


    @unittest.skip("Disable if too heavy")
    def test_c_and_sklearn_over_and_over(self):
        np.random.seed(0)
        for i in range(100):
            self.test_c_and_sklearn_equal()
            #self.test_py_and_mine_equal()
            pass


    @unittest.skip("enable later")
    def test_py_and_c_equal_random(self):
        print("test_py_and_c_equal_random() - start")
        fit_params = list(randomize_fit_params())
        centroids_list_py = kmeans_1.KmeansAlgorithm(*fit_params)
        centroids_list_c = mykmeanssp.fit(*fit_params)
        centroids_list_py = np.array(centroids_list_py)
        centroids_list_c = np.array(centroids_list_c)
        dist_c  = np.all(np.abs(centroids_list_py-centroids_list_c) < 0.001)
        assert(dist_c)
    
    @unittest.skip("Passing")
    def check_crash(self):
        with open('/home/ubuntu/lala2.bin','rb') as f:
            fit_params = pickle.loads(f.read())
        centroids_list_py = kmeans_pp.KmeansAlgorithm(*fit_params)
        print("survived")
    
    #@unittest.skip("What do we need this for")
    def test_my_py_runtime_vs_sklearn(self):
        print("Blargh")
        for i in range(1000):
            np.random.seed(i)
            #print(f"hai {i}")
            fit_params = list(randomize_fit_params(k=10, max_iter=1000, eps=None, point_count=150, dims_count=7))
            with open('/home/ubuntu/lala2.bin','wb') as f:
                f.write(pickle.dumps(fit_params))
            time_1, time_2, time_3 = time.time(), time.time(), time.time()
            time_1 = time.time()
            centroids_list_py = kmeans_pp.KmeansAlgorithm(*fit_params)
            time_2 = time.time()
            centroids_list_sk = kmeans_sk.KmeansAlgorithm(*fit_params)
            time_3 = time.time()
            delta_py = time_2-time_1
            delta_sk = time_3-time_2
            relative_sk_py = relative_error_centroids(centroids_list_sk, centroids_list_py)
            #print(f"delta_py/delta_sk={delta_py/delta_sk}")
            print(f"relative err = {relative_sk_py}")
            pass

    @unittest.skip("Only needed once in a while")
    def test_equal_to_templates(self):
        def test_equal_to_template_idx(*args):
            print(f"test_equal_to_template_idx{args} - start")
            with patch('sys.argv', ["python3","blah.py"]+list([str(x) for x in args])):
                fit_params = list(kmeans_pp.extract_fit_params())
                initial_centroids_list = fit_params[0]

                centroids_list_desired = None
                file_expected = args[-1].replace("input","output")
                file_expected = re.sub('_db_\d', '', file_expected)
                with open(file_expected, 'r') as f:
                    s = f.read().split("\n")[:-1]
                    s = [x.split(",") for x in s]
                    s[0] = [int(y) for y in s[0]]
                    initial_centroids_list_desired = s[0]
                    centroids_list_desired = np.array([[float(y) for y in x] for x in s[1:]])

                self.assertTrue(np.all(np.array(initial_centroids_list)==np.array(initial_centroids_list_desired)))

                centroids_list_py = kmeans_pp.KmeansAlgorithm(*fit_params)
                centroids_list_sk = kmeans_sk.KmeansAlgorithm(*fit_params)

                centroids_list_py = np.array(centroids_list_py)
                centroids_list_sk = np.array(centroids_list_sk)

                dist_desired_py = dist_between_centroid_lists(centroids_list_desired, centroids_list_py)
                dist_py_sk      = dist_between_centroid_lists(centroids_list_py, centroids_list_sk)
                
                print(dist_desired_py)
                print(dist_py_sk)

                self.assertTrue(dist_desired_py == 0)
                #self.assertTrue(dist_desired_sk == 0)
    
        print("test_equal_to_templates() - start")
        pathloc = f"/home/ubuntu/repos/softproj_2/resources/test_data_2"
        test_equal_to_template_idx(3,  333,                       0, f"{pathloc}/input_1_db_1.txt", f"{pathloc}/input_1_db_2.txt")
        test_equal_to_template_idx(7,  kmeans_pp.MAX_ITER_UNSPEC, 0, f"{pathloc}/input_2_db_1.txt", f"{pathloc}/input_2_db_2.txt")
        test_equal_to_template_idx(15, 750,                       0, f"{pathloc}/input_3_db_1.txt", f"{pathloc}/input_3_db_2.txt")


def randomize_fit_params(k=None, max_iter=None, eps=None, point_count=None, dims_count=None):
    #np.random.seed(4)
    k = k or np.random.randint(2, 5)
    max_iter = max_iter or np.random.randint(100, 300)
    eps = eps or float(np.random.rand(1,1))/100
    point_count = point_count or np.random.randint(50, 200)
    dims_count = dims_count or np.random.randint(3,5)
    datapoints_list = np.random.rand(point_count, dims_count)
    indices_column = np.array(range(datapoints_list.shape[0])).reshape(datapoints_list.shape[0],1)
    datapoints_list = np.hstack((indices_column, datapoints_list))
    initial_centroids_list = kmeans_pp.KMeansPlusPlus(k, datapoints_list)
    initial_centroids_indices_as_written = [int(initial_centroids_list[i][0]) for i in range(len(initial_centroids_list))]
    #print(','.join([str(x) for x in initial_centroids_indices_as_written]))
    initial_centroids_indices_actual = kmeans_pp.select_actual_centroids(datapoints_list, initial_centroids_list)
    datapoints_list = datapoints_list[:,1:]
    return (
        initial_centroids_indices_actual,
        datapoints_list,
        dims_count,
        k,
        point_count,
        max_iter,
        eps
    )


def idx_of_centroid_closest_to_point(centroids_list: List[List[float]], point: List[float]):
    dist_vector = np.sum(np.square(np.array(centroids_list)-np.array(point)), axis=1)
    return (np.argmin(dist_vector), np.min(dist_vector))

def dist_between_centroid_lists_redundant(list_1: List[List[float]], list_2: List[List[float]]) -> float:
    if len(list_1) != len(list_2): return math.inf
    k = len(list_1)
    if k == 0: return 0
    dims = len(list_1[0])
    for i in range(k):
        if (len(list_1[i]) != dims) or (len(list_2[i]) != dims): return math.inf
    dist = 0
    for i in range(k):
        #print(f"UNO: list_1 length is {list_1}, list_2 is {(list_2)}")
        idx, _dist = idx_of_centroid_closest_to_point(list_1, list_2[0])
        dist += _dist
        list_1 = np.array(list(list_1)[:idx] + list(list_1)[idx+1:])
        list_2 = np.array(list(list_2)[1:])
        #print(f"list_1 length is {list_1}, list_2 is {(list_2)}")
    return dist

def dist_between_centroid_lists(list_1: np.ndarray, list_2: np.ndarray) -> float:
    return np.linalg.norm((list_1)-(list_2))

def relative_error_centroids(centroids_real: List[List[float]], centroids_calc: List[List[float]]) -> float:
    return np.linalg.norm(np.sort(centroids_real)-np.sort(centroids_calc))/np.linalg.norm(centroids_real)

if __name__ == '__main__':
    print("Starting tests")
    unittest.main()
