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

print(result)
