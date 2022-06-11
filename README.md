This project is an assignment I was given in one of my courses

It performs Kmeans clustering on a given dataset (csv/txt) using the Kmeans++ initialization method, using Numpy, CPython, and testing against Sklearn

I implemented the algorithm in bare C, and wrote a wrapper in CPython for easy interfacing with Python - see `source/kmeans.c`

The initialization (Kmeans++) is performed in `source/kmeans_pp.py` (which is also where you'll find the main thread of the program) in `KMeansPlusPlus()`.

For the initialization part, there is also a more elegant solution in `KMeansPlusPlus_original()`, however it randomizes in a way that is slightly different than the one provided in the example outputs, so I ended up not using it.

I also wrote some pretty extensive tests - among other things, comparing the results to sklearn, in `source/tests.py`

* I am aware of the mess in the file structure, also I would usually try to modularize the files more, but the project required the files to have specific names

* I also didn't pick the uninformative error messages
