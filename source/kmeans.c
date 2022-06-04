#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h> 


/*/////////////////////
// Constants - START //
/////////////////////*/

#define MAX_ITER_UNSPEC 200

#define MSG_ERR_INVALID_INPUT "Invalid Input!\n"
#define MSG_ERR_GENERIC       "An Error Has Occurred\n"

#define STATUS_SUCCESS      0
#define STATUS_ERROR        1
#define STATUS_ERROR_MALLOC STATUS_ERROR
#define STATUS_ERROR_FOPEN  STATUS_ERROR
#define STATUS_ERROR_FORMAT 2

#define RESULT_FOPEN_SUCCESS 0
#define RESULT_FOPEN_ERROR   1

#define EPSILON ((double)0.001)

/*xxxxxxxxxxxxxxxxxxx//
//  Constants - END  //
//xxxxxxxxxxxxxxxxxxx*/


/*///////////////////////////
// Data Structures - START //
///////////////////////////*/

typedef struct Point {
    struct Point* next;
    double* coord;
    int cluster;
} point_t;

typedef struct PointList {
    point_t* head;
    point_t* tail;
    int count;
    int dims;
} pointlist_t;

point_t* point_init(int dims) {
    point_t* point;
    point = malloc(sizeof(point_t));
    if (!point) return NULL;
    point->next = NULL;
    point->coord = malloc(sizeof(double)*dims);
    if (!point->coord) {
        free(point);
        return NULL;
    }
    point->cluster = -1;
    return point;
}

void point_free(point_t** point, int recursive) {
    point_t* next;
    next = (*point)->next;
    if (!(*point)) return;
    if ((*point)->coord) free((*point)->coord);
    free((*point));
    if (recursive) point_free(&next, 1);
}

pointlist_t* pointlist_init(int dims) {
    pointlist_t* list;
    list = malloc(sizeof(pointlist_t));
    if (!list) return NULL;
    list->head = NULL;
    list->tail = NULL;
    list->dims = dims;
    list->count = 0;
    return list;
}

/*void pointlist_free(pointlist_t** list) {
    if (!(*list)) return;
    if ((*list)->head) {
        point_free(&((*list)->head), 1);
    }
    free((*list));
}*/

void pointlist_free(point_t** list, int count, int dims) {
    int i;
    for (i=0; i<count; i++) {
        if((*list)[i].coord) free((*list)[i].coord);
    }
    free((*list));
}

/*xxxxxxxxxxxxxxxxxxxxxxxxx//
//  Data Structures - END  //
//xxxxxxxxxxxxxxxxxxxxxxxxx*/


/*///////////////////////////////
// Method Declarations - START //
///////////////////////////////*/

int allocate_centroids(point_t** centroids_list, int k, int dims_count);
void init_centroids(point_t* centroids_list, point_t* points_list, int k, int dims_count);
void assign_every_point_to_nearest_cluster(point_t* centroids_list, point_t* points_list, int k, int line_count, int dims_count);
void copy_centroids(point_t* dst, point_t* src, int k, int dims_count);
void reassign_all_centroids(point_t* centroids_list, point_t* points_list, int k, int line_count, int dims_count);
void kmeans_iteration(point_t* centroids_list, point_t* points_list, int k, int line_count, int dims_count);
int is_convergence(point_t* prev_centroids_list, point_t* centroids_list, int k, int dims_count);
int write_results_to_file(char* path_to_output, point_t* centroids_list, int k, int dims_count);

/* point_t methods */
double get_dist_between_points(point_t point_1, point_t point_2, int dims_count);
double get_abs_point(point_t point, int dims_count);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//  Method Declarations - END  //
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/*//////////////////////////////
// Method Definitions - START //
//////////////////////////////*/

double get_square_dist_between_points(point_t point_1, point_t point_2, int dims_count) {
    int d;
    double dist;
    dist = 0;
    for (d = 0; d < dims_count; d++) {
        dist += pow((point_1.coord[d] - point_2.coord[d]), 2);
    }
    return dist;
}

double get_dist_between_points(point_t point_1, point_t point_2, int dims_count) {
    return pow(get_square_dist_between_points(point_1, point_2, dims_count), 0.5);
}

double get_abs_point(point_t point, int dims_count) {
    int d;
    double dist;
    dist = 0;
    for (d = 0; d < dims_count; d++) {
        dist += pow(point.coord[d], 2);
    }
    dist = pow(dist, 0.5);
    return dist;
}

int isanum(char* s) {
    int i;
    for (i = 0; s[i] != 0; i++) {
        if (!(('0' <= s[i]) && (s[i] <= '9')))
            return 0;
    }
    return 1;
}

int allocate_points(point_t** centroids_list, int points_count, int dims_count) {
    int i, j;
    (*centroids_list) = malloc(points_count * sizeof(point_t));
    if ((*centroids_list) == NULL) {
        return STATUS_ERROR_MALLOC;
    }
    for (i = 0; i < points_count; i++) {
        (*((*centroids_list) + i)).coord   = malloc(dims_count * sizeof(double));
        (*((*centroids_list) + i)).cluster = -1;
        if ((*((*centroids_list) + i)).coord == NULL) {
            for (j = 0; j < i; j++) {
                free((*centroids_list)[j].coord);
            }
            free(*centroids_list);
            return STATUS_ERROR_MALLOC;
        }
    }
    return STATUS_SUCCESS;
}

int allocate_centroids(point_t** centroids_list, int k, int dims_count) {
    int i, status;
    status = allocate_points(centroids_list, k, dims_count);
    if (status != STATUS_SUCCESS) return status;
    for (i=0; i<k; i++) {
        (*centroids_list)[i].cluster = i;
    }
    return STATUS_SUCCESS;
}

void init_centroids(point_t* centroids_list, point_t* points_list, int k, int dims_count) {
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < dims_count; j++) {
            centroids_list[i].coord[j] = points_list[i].coord[j];
        }
        points_list[i].cluster = i;
    }
}

void assign_every_point_to_nearest_cluster(point_t* centroids_list, point_t* points_list, int k, int line_count, int dims_count) {
    int i, j, min_j;
    double dist, min_dist;
    for (i = 0; i < line_count; i++) {
        min_j = 0;
        min_dist = get_square_dist_between_points(points_list[i], centroids_list[0], dims_count);
        for (j = 1; j < k; j++) {
            dist = get_square_dist_between_points(points_list[i], centroids_list[j], dims_count);
            if (dist < min_dist) {
                min_j = j;
                min_dist = dist;
            }
        }
        points_list[i].cluster = min_j;
    }
}

void copy_centroids(point_t* dst, point_t* src, int k, int dims_count) {
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < dims_count; j++) {
            dst[i].coord[j] = src[i].coord[j];
        }
    }
}

void reassign_all_centroids(point_t* centroids_list, point_t* points_list, int k, int line_count, int dims_count) {
    int i, j, d;
    int number_of_hits;
    number_of_hits = 0;
    for (i = 0; i < k; i++) {
        number_of_hits = 0;
        for (j = 0; j < k; j++) centroids_list[i].coord[j] = 0;
        for (j = 0; j < line_count; j++) {
            if (points_list[j].cluster == i) {
                number_of_hits++;
                for (d = 0; d < dims_count; d++)
                    centroids_list[i].coord[d] += points_list[j].coord[d];
            }
        }
        if (number_of_hits > 0) {
            for (j = 0; j < dims_count; j++) {
                centroids_list[i].coord[j] /= number_of_hits;
            }
        }
    }
}

void kmeans_iteration(point_t* centroids_list, point_t* points_list, int k, int line_count, int dims_count) {
    assign_every_point_to_nearest_cluster (centroids_list, points_list, k, line_count, dims_count);
    reassign_all_centroids                (centroids_list, points_list, k, line_count, dims_count);
}

int is_convergence(point_t* prev_centroids_list, point_t* centroids_list, int k, int dims_count) {
    int i;
    double dist;
    for (i = 0; i < k; i++) {
        dist = get_dist_between_points(prev_centroids_list[i], centroids_list[i], dims_count);
        if (dist >= EPSILON) return 0;
    }
    return 1;
}

int write_results_to_file(char* path_to_output, point_t* centroids_list, int k, int dims_count) {
    FILE* fh;
    int i, j;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    errno_t err;
    err = fopen_s(&fh, path_to_output, "wb");
    if (err != RESULT_FOPEN_SUCCESS) {
        return STATUS_ERROR_FOPEN;
    }
#else
    fh = fopen(path_to_output, "wb");
    if (fh == NULL) {
        return STATUS_ERROR_FOPEN;
    }
#endif


    for (i = 0; i < k; i++) {
        for (j = 0; j < dims_count-1; j++)
            fprintf(fh, "%.4f,", centroids_list[i].coord[j]);
        fprintf(fh, "%.4f\n", centroids_list[i].coord[dims_count - 1]);
    }

    fclose(fh);
    return STATUS_SUCCESS;
}

void print_centroids(point_t* centroids_list, int k, int dims_count) {
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < dims_count - 1; j++)
            printf("%.4f,", centroids_list[i].coord[j]);
        printf("%.4f\n", centroids_list[i].coord[dims_count - 1]);
    }
    printf("\n");
}

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//  Method Definitions - END  //
//xxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


static int set_points_to_dst(PyObject* obj, int point_count, int dims_count, point_t* dst) {
    PyObject *points_iter, *next_obj, *coords_iter, *next_coord;
    int i, j;
    
    /*
    *dst = malloc(sizeof(point_t)*point_count);
    for (i=0; i<point_count; i++) {
        (*dst)[i].coords = malloc(sizeof(double)*dims_count);
        if (!(*dst)[i].coords) {
            for (j=0; j<i; j++) free((*dst)[j].coords);
            free(*dst);
            *dst = NULL;
        }
    }*/

    points_iter = PyObject_GetIter(obj);

    i=0;
    while ((next_obj = PyIter_Next(points_iter))) {

        coords_iter = PyObject_GetIter(next_obj);
        j = 0;
        while ((next_coord = PyIter_Next(coords_iter))) {
            (/* * */dst)[i].coord[j] = PyFloat_AsDouble(next_coord);
            j++;
            Py_DECREF(next_coord);
        }
        Py_DECREF(coords_iter);
        Py_DECREF(next_obj);
        i++;
    }
    Py_DECREF(points_iter);
    
    return 0;
}

static void get_points_list(PyObject* obj_initial_centroids, PyObject* obj_datapoints,
                                int initial_centroids_count, int datapoints_count,
                                int dims_count,
                                point_t* initial_centroids, point_t* datapoints) {
    int i;
    set_points_to_dst(obj_initial_centroids, initial_centroids_count, dims_count, initial_centroids);
    set_points_to_dst(obj_datapoints, datapoints_count, dims_count, datapoints);
    for (i=0; i<initial_centroids_count; i++) {
        (/* * */initial_centroids)[i].cluster = i;
    }
}

static int allocate_necessary(point_t** centroids_list, point_t** prev_centroids_list, point_t** points_list,
                            int dims_count, int k, int point_count) {
    int status;

    status = allocate_points(points_list, point_count, dims_count);
    if (status != STATUS_SUCCESS) {
        return status;
    }

    /* allocate memory for centroids */
    status = allocate_centroids(centroids_list, k, dims_count);
    if (status != STATUS_SUCCESS) {
        pointlist_free(points_list, point_count, dims_count);
        return status;
    }

    /* allocate memory for old copy of centroids (needed to calculate convergence) */
    status = allocate_centroids(prev_centroids_list, k, dims_count);
    if (status != STATUS_SUCCESS) {
        pointlist_free(points_list, point_count, dims_count);
        pointlist_free(centroids_list, k, dims_count);
        return status;
    }

    return STATUS_SUCCESS;
}

static void reach_convergence(point_t* centroids_list, point_t* prev_centroids_list, point_t* points_list,
                            int dims_count, int k, int point_count, int max_iter) {
    int i;
    max_iter = 0; /* FIXME */
    /* perform kmeans algorithm! */
    /* init_centroids(centroids_list, points_list, k, dims_count); */
    for (i = 0; i < max_iter; i++) {
        copy_centroids(prev_centroids_list, centroids_list, k, dims_count);
        kmeans_iteration(centroids_list, points_list, k, point_count, dims_count);
        if (is_convergence(prev_centroids_list, centroids_list, k, dims_count)) break;
    }
}

static point_t* calculate_centroids(PyObject* obj_initial_centroids, PyObject* obj_datapoints,
                                    int dims_count, int k, int point_count, int max_iter) {
    int status;
    point_t *centroids_list, *points_list, *prev_centroids_list;

    status = allocate_necessary(&centroids_list, &prev_centroids_list, &points_list,
                        dims_count, k, point_count);
    if (status != STATUS_SUCCESS) {
        return NULL;
    }

    /* assign points from py objs */
    get_points_list(obj_initial_centroids, obj_datapoints,
                    k, point_count,
                    dims_count,
                    centroids_list, points_list);
    
    reach_convergence(centroids_list, prev_centroids_list, points_list,
                        dims_count, k, point_count, max_iter);
    
    pointlist_free(&points_list, point_count, dims_count);
    pointlist_free(&prev_centroids_list, k, dims_count);

    return centroids_list;
}

static PyObject* centroids_to_PyObject(point_t* centroids_list, int k, int dims_count) {
    double z = 1;
    int n = 2;

    printf("Got length %d, dims %d\n", k, dims_count);

    return Py_BuildValue("d", 5.3);
}

static PyObject* fit(PyObject *self, PyObject *args) {
    point_t *centroids_list;
    PyObject *obj_initial_centroids, *obj_datapoints;
    PyObject *centroids_list_as_pyobject;
    int point_count, dims_count;
    int k;
    int max_iter;

    if (!PyArg_ParseTuple(args, "OOiiii", &obj_initial_centroids, &obj_datapoints, &dims_count, &k, &point_count, &max_iter)) {
        return NULL;
    }

    centroids_list = calculate_centroids(obj_initial_centroids, obj_datapoints, dims_count, k, point_count, max_iter);
    centroids_list_as_pyobject = centroids_to_PyObject(centroids_list, k, dims_count);
    pointlist_free(&centroids_list, k, dims_count);

    return centroids_list_as_pyobject;

}

static PyMethodDef mykmeansspMethods[] = {
    {"fit", 
      (PyCFunction) fit,
      METH_VARARGS, 
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")},
    {NULL, NULL, 0, NULL} 
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    mykmeansspMethods /* the PyMethodDef array from before containing the methods of the extension */
};


/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
