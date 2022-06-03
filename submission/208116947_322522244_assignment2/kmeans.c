/*
// kmeans.c
// scroll down all the way to the end to see main()
*/

/*////////////////////
// Includes - START //
////////////////////*/

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*xxxxxxxxxxxxxxxxxx//
//  Includes - END  //
//xxxxxxxxxxxxxxxxxx*/


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
    double* coord;
    int cluster;
} point_t;

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

/* File logistics */
int FILE_locate(FILE* fh, char needle);
char* FILE_get(FILE* fh, int length);
char* FILE_get_next_line(FILE* fh);
char* FILE_get_next_num(FILE* fh);
int get_number_of_dimensions(FILE* fh);
int get_number_of_lines(FILE* fh);

/*
// INPUT  : argc, argv
// OUTPUT : k, max_iter, path_to_input, path_to_output
// RETVAL : STATUS_SUCCESS or STATUS_ERROR
*/
int set_arguments(int argc, char* argv[],
        int* k,
        int* max_iter,
        char** path_to_input,
        char** path_to_output);

/*
// INPUT  : path_to_input
// OUTPUT : k, max_iter, path_to_input, path_to_output
// RETVAL : STATUS_SUCCESS or STATUS_ERROR
*/
int read_data(char* path_to_input,
        point_t** points_list, int* line_count, int* dims_count);

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

int set_arguments(int argc, char* argv[],
        int* k,
        int* max_iter,
        char** path_to_input,
        char** path_to_output) {

    char* tmp__string_part;

    switch (argc) {
        case 4: {
            if (!isanum(argv[1])) return STATUS_ERROR;
            *k = (int)strtol(argv[1], &tmp__string_part, 10);
            if ((*k) < 1) return STATUS_ERROR;
            *max_iter = MAX_ITER_UNSPEC;
            *path_to_input = argv[2];
            *path_to_output = argv[3];
            return STATUS_SUCCESS;
        }
        case 5: {
            if (!isanum(argv[1]) || !isanum(argv[2])) return STATUS_ERROR;
            *k = (int)strtol(argv[1], &tmp__string_part, 10);
            if ((*k) < 1) return STATUS_ERROR;
            *max_iter = (int)strtol(argv[2], &tmp__string_part, 10);
            *path_to_input = argv[3];
            *path_to_output = argv[4];
            return STATUS_SUCCESS;
        }
        default: {
            return STATUS_ERROR;
        }
    }
}

int FILE_locate(FILE* fh, char needle) {
    char c;
    int counted;
    long start;
    start = ftell(fh);
    counted = 0;
    while (1) {
        c = fgetc(fh);
        if ((c == EOF) || (c == needle)) {
            break;
        }
        counted++;
    }
    fseek(fh, start, SEEK_SET);

    if (c != needle) return -1;
    else return counted;
}

char* FILE_get(FILE* fh, int length) {
    char* line;
    if (length < 0) length = 0;
    line = malloc(/*sizeof(char) * */(length + 1));
    if (line == NULL) return NULL;

    fgets(line, length+1, fh);

    return line;
}

char* FILE_get_next_line(FILE* fh) {
    char* line;
    int length = FILE_locate(fh, '\n');
    if (length < 0) length = 0;
    line = FILE_get(fh, length);
    if ((line != NULL) && (!feof(fh))) fgetc(fh); /* drop last */
    return line;
}

char* FILE_get_next_num(FILE* fh) {
    int next_comma, next_newl;
    char* line;

    next_comma = FILE_locate(fh, ',');
    next_newl  = FILE_locate(fh, '\n');
    if ((next_comma == -1) || ((next_newl != -1) && (next_newl < next_comma))) {
        next_comma = next_newl;
    }

    if (next_comma < 0) next_comma = 0;


    line = FILE_get(fh, next_comma);

    if ((line != NULL) && (!feof(fh))) fgetc(fh); /* drop last */

    return line;
}

int get_number_of_dimensions(FILE* fh) {
    int i, d, got_dim;
    char* line;
    d = 1;
    got_dim = 0;
    while (!feof(fh)) {
        line = FILE_get_next_line(fh);
        if ((line[0] == 0) || (line[0] == '\n') || (line[0] == ',')) {
            free(line);
            continue;
        }
        got_dim = 1;
        for (i = 0; (line[i] != 0) && (line[i] != '\n'); i++) {
            if (line[i] == ',') d++;
        }
        free(line);
        break;
    }
    if (!got_dim) d = 0;
    fseek(fh, 0, SEEK_SET);
    return d;
}

int get_number_of_lines(FILE* fh) {
    int lines;
    char c;
    lines = 0;
    while (1) {
        c = fgetc(fh);
        if (c == EOF) break;
        if (c == '\n') lines++;
    }
    fseek(fh, 0, SEEK_SET);
    return lines;
}

int read_data(char* path_to_input,
        point_t** points_list, int* line_count, int* dims_count) {
    int i, j;
    char* line;
    FILE* fh;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    errno_t err;
    err = fopen_s(&fh, path_to_input, "rb");
    if (err != RESULT_FOPEN_SUCCESS) {
        return STATUS_ERROR_FOPEN;
    }
#else
    fh = fopen(path_to_input, "rb");
    if (fh == NULL) {
        return STATUS_ERROR_FOPEN;
    }
#endif

    /* we assume that every line gets the same number of commas */
    *dims_count = get_number_of_dimensions(fh);
    *line_count = get_number_of_lines(fh);
    if ((*dims_count <= 0) || (*line_count <= 0)) {
        fclose(fh);
        return STATUS_ERROR_FORMAT;
    }

    (*points_list) = malloc(*line_count * sizeof(point_t));
    if ((*points_list) == NULL) return STATUS_ERROR_MALLOC;
    for (i = 0; i < *line_count; i++) {
        (*points_list)[i].coord   = malloc(*dims_count * sizeof(double));
        (*points_list)[i].cluster = -1;
        if ((*points_list)[i].coord == NULL) {
            for (j = 0; j < i; j++) {
                free((*points_list)[j].coord);
            }
            free(*points_list);
            return STATUS_ERROR_MALLOC;
        }

        for (j = 0; j < *dims_count; j++) {
            line = FILE_get_next_num(fh);
            ((*points_list)[i]).coord[j] = atof(line);
            free(line);
        }
    }

    fclose(fh);
    if ((*points_list) == NULL) return STATUS_ERROR_MALLOC;
    return STATUS_SUCCESS;
}

int allocate_centroids(point_t** centroids_list, int k, int dims_count) {
    int i, j;
    (*centroids_list) = malloc(k * sizeof(point_t));
    if ((*centroids_list) == NULL) {
        return STATUS_ERROR_MALLOC;
    }
    for (i = 0; i < k; i++) {
        (*((*centroids_list) + i)).coord   = malloc(dims_count * sizeof(double));
        (*((*centroids_list) + i)).cluster = i;
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


int main(int argc, char* argv[]) {
    
    /* input arguments */
    int k;
    int max_iter;
    char* path_to_input;
    char* path_to_output;
    point_t* points_list;
    point_t* centroids_list;
    point_t* prev_centroids_list;
    /* results from read_data */
    int line_count, dims_count;
    /* misc */
    int error;
    int status;
    int i;

    error = 0; /* innocent until proven guilty */

    /* set input arguments, or die if error */
    if (STATUS_SUCCESS != 
            set_arguments(
                argc, argv,
                &k, &max_iter, &path_to_input, &path_to_output)) {
        printf(MSG_ERR_INVALID_INPUT);
        return 1;
    }

    /* read input file to points_list */
    status = read_data(path_to_input, &points_list, &line_count, &dims_count);
    if (status != STATUS_SUCCESS) {
        error = 1;
        goto right_before_return;
    }

    if (line_count < k) {
        error = 1;
        goto free_points;
    }

    /* allocate memory for centroids */
    status = allocate_centroids(&centroids_list, k, dims_count);
    if (status != STATUS_SUCCESS) {
        error = 1;
        goto free_points;
    }

    /* allocate memory for old copy of centroids (needed to calculate convergence) */
    status = allocate_centroids(&prev_centroids_list, k, dims_count);
    if (status != STATUS_SUCCESS) {
        error = 1;
        goto free_centroids;
    }


    /* perform kmeans algorithm! */
    init_centroids(centroids_list, points_list, k, dims_count);
    for (i = 0; i < max_iter; i++) {
        copy_centroids(prev_centroids_list, centroids_list, k, dims_count);
        kmeans_iteration(centroids_list, points_list, k, line_count, dims_count);
        if (is_convergence(prev_centroids_list, centroids_list, k, dims_count)) break;
    }

    /* once converged or reached max iterations, attempt to print to file */
    status = write_results_to_file(path_to_output, centroids_list, k, dims_count);
    if (status != STATUS_SUCCESS) {
        error = 1;
        goto free_prev_centroids;
    }


    /* cleanup */

    free_prev_centroids:
    for (i = 0; i < k; i++) {
        free((prev_centroids_list[i]).coord);
    }
    free(prev_centroids_list);

    free_centroids:
    for (i = 0; i < k; i++) {
        free((centroids_list[i]).coord);
    }
    free(centroids_list);

    free_points:
    for (i = 0; i < line_count; i++) {
        free((points_list[i]).coord);
    }
    free(points_list);

    right_before_return:
    if (error == 1) printf(MSG_ERR_GENERIC);
    return error;
}
