#!/usr/bin/env python3

import unittest
from unittest.mock import patch
import logging
import sys
import numpy as np
import kmeans_pp
import kmeans as kmeans_1
import mykmeanssp


class TestFit(unittest.TestCase):

    @unittest.skip("Only needed once in a while")
    def test_equal_to_templates(self):
        def test_equal_to_template_idx(*args):
            print("test_equal_to_template_idx() - start")
            with patch('sys.argv', ["python3","blah.py"]+list([str(x) for x in args])):
                print(f"Sys argv: {sys.argv}")
                fit_params = list(kmeans_pp.extract_fit_params())
                desired = None
                file_expected = args[-1].split("_")[0].replace("input","output")+"_"+args[-1].split("_")[1]+".txt"
                with open(file_expected, 'r') as f:
                    s = f.read().split("\n")[:-1]
                    s = [x.split(",") for x in s]
                    s[0] = [int(y) for y in s[0]]
                    initial_centroids_list = s[0]
                    desired = np.sort(np.array([[float(y) for y in x] for x in s[1:]]))
                fit_params[0] = initial_centroids_list
                centroids_list_py = kmeans_1.KmeansAlgorithm(*fit_params)
                centroids_list_c = mykmeanssp.fit(*fit_params)
                centroids_list_py = np.sort(np.array(centroids_list_py))
                centroids_list_c = np.sort(np.array(centroids_list_c))
                dist_py = np.all(np.abs(desired-centroids_list_py) < 0.001)
                assert(dist_py)
                dist_c  = np.all(np.abs(centroids_list_py-centroids_list_c) < 0.001)
                assert(dist_c)
        print("test_equal_to_templates() - start")
        test_equal_to_template_idx(3, 333, 0, "input_1_db_1.txt", "input_1_db_2.txt")
        test_equal_to_template_idx(7, kmeans_pp.MAX_ITER_UNSPEC, 0, "input_2_db_1.txt", "input_2_db_2.txt")
        test_equal_to_template_idx(15, 750, 0, "input_3_db_1.txt", "input_3_db_2.txt")


    def test_py_and_c_equal_random(self):
        print("test_py_and_c_equal_random() - start")
        def randomize_fit_params():
            k, max_iter, eps = np.random.randint(2, 10), np.random.randint(100, 300), float(np.random.rand(1,1))/100
            point_count = np.random.randint(50, 200)
            dims_count = np.random.randint(3,6)
            datapoints_list = list(np.random.rand(point_count, dims_count))
            a = np.array(range(point_count))[...,None]
            datapoints_list = np.hstack((a, datapoints_list))
            initial_centroids_list = kmeans_pp.KmeansPlusPlus(k, datapoints_list)
            datapoints_list = sorted(datapoints_list, key=lambda x: float(x[0]))
            datapoints_list = [x[1:] for x in datapoints_list]

            return (
                initial_centroids_list,
                datapoints_list,
                dims_count,
                k,
                point_count,
                max_iter,
                eps
            )
        fit_params = list(randomize_fit_params())
        centroids_list_py = kmeans_1.KmeansAlgorithm(*fit_params)
        centroids_list_c = mykmeanssp.fit(*fit_params)
        centroids_list_py = np.sort(np.array(centroids_list_py))
        centroids_list_c = np.sort(np.array(centroids_list_c))
        dist_c  = np.all(np.abs(centroids_list_py-centroids_list_c) < 0.001)
        assert(dist_c)


if __name__ == '__main__':
    print("Starting tests")
    unittest.main()