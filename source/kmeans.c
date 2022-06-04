#define PY_SSIZE_T_CLEAN
#include <Python.h> 
#include <math.h> 


static double geo_c(double z, int n)
{
    double geo_sum = 0;
    int i;
    for (i=0; i<n; i++){
        /* pow(x,y) function raises x to the power of y - it is from <math.h> */
        geo_sum += pow(z,i);
     }
    printf("Blargh\n");
    return geo_sum;
}


static PyObject* geo_capi(PyObject *self, PyObject *args)
{
    double z;
    int n;
    if(!PyArg_ParseTuple(args, "di", &z, &n)) {
        return NULL;
    }

    return Py_BuildValue("d", geo_c(z, n));
}

static PyMethodDef capiMethods[] = {
    {"geo", 
      (PyCFunction) geo_capi,
      METH_VARARGS, 
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")},
    {NULL, NULL, 0, NULL} 
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "capi_demo1", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};


/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
