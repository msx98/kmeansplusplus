#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void pointlist_free(pointlist_t** list) {
    if (!(*list)) return;
    if ((*list)->head) {
        point_free(&((*list)->head), 1);
    }
    free((*list));
}

point_t* pointlist_add(pointlist_t* list, double* coord) {
    int i;
    list->tail->next = point_init(list->dims);
    if (!list->tail->next) return NULL;
    for (i=0; i<list->dims; i++) {
        list->tail->next->coord[i] = coord[i];
    }
}


typedef struct Node {
    struct Node* next;
    void* data;
} node_t;

typedef struct List {
    node_t* first;
    node_t* last;
    int length;
    size_t sizeof_node;
    int contains_lists;
} list_t;

node_t* node_init(void* data, size_t sizeof_data) {
    int i;
    node_t* head;
    head = malloc(sizeof(node_t));
    if (!head) return NULL;
    head->data = NULL;
    head->data = malloc(sizeof_data);
    if (!head->data) {
        free(head);
        return NULL;
    }
    for (i=0; i<sizeof_data; i++) {
        *((char *)(head->data  + i)) =  *((char *)(data + i));
    }
    head->next = NULL;
    return head;
}

void node_free(node_t** node) {
    if ((*node)->data) {
        free((*node)->data);
    }
    free(*node);
}

list_t* list_init(size_t sizeof_node) {
    list_t* list;
    list = malloc(sizeof(list_t));
    if (!list) return NULL;
    list->first = NULL;
    list->last = NULL;
    list->length = 0;
    if (sizeof_node == 0) {
        list->sizeof_node = sizeof(list_t);
        list->contains_lists = 1;
    } else {
        list->sizeof_node = sizeof_node;
        list->contains_lists = 0;
    }
    return list;
}

node_t* list_append(list_t* list, void *data) {
    node_t *node, *last;
    node = node_init(data, list->sizeof_node);
    if (!node) return NULL;
    if (list->length == 0) {
        list->first = node;
        list->last  = node;
    } else {
        list->last->next = node;
        list->last = node;
    }
    list->length++;
    return node;
}

node_t* list_get(list_t* list, int idx) {
    node_t* p = list->first;
    while (p && (idx > 0)) {
        p = p->next;
        idx--;
    }
    return p;
}

node_t* list_pop(list_t* list, int idx) {
    node_t *node_before, *node, *node_after;
    if ((idx < 0) || (idx >= list->length)) {
        /* index is invalid */
        return NULL;
    }
    
    /* index is valid */

    if (list->length == 1) {
        node = list->first;
        list->first = NULL;
        list->last = NULL;
        return node;
    }

    /* index is valid, list length >= 2 */

    if (idx == 0) {
        node = list->first;
        list->first = list->first->next;
        return node;
    }

    /* index >= 1 is valid, list length >= 2 */

    node_before = list_get(list, idx-1);
    if (idx == list->length-1) {
        list->last = node_before;
    }
    node = node_before->next;
    node_before->next = NULL;
    list->length--;
    return node;
}

void list_free(list_t** list) {
    node_t *p, *p_next;
    list_t *inner_list;
    p = (*list)->first;
    
    if ((*list)->contains_lists) {
        while (p) {
            p_next = p->next;
            inner_list = (list_t*) p->data;
            if (inner_list) list_free(&inner_list);
            free(p);
            p = p_next;
        }
    } else {
        while (p) {
            p_next = p->next;
            if (p->data) free(p->data);
            free(p);
            p = p_next;
        }
    }
    free(*list);
}

node_t* push2(node_t* last, void *x, size_t x_size){
    
    node_t* new_node = (node_t*)malloc(sizeof(node_t));
 
    new_node->data  = malloc(x_size);
    if (!new_node->data) {
        free(new_node);
        return NULL;
    }

    new_node->next = NULL;
        
    // here we are copying the data x to new node created
    for(int i=0; i < x_size; i++){
        *(char *)(new_node->data + i) = *(char *)(x + i);
    }

    last->next = new_node;

    return new_node;
}

node_t* node2_init(void *x, size_t x_size) {
    node_t* new_node = (node_t*)malloc(sizeof(node_t));
 
    new_node->data  = malloc(x_size);
    if (!new_node->data) {
        free(new_node);
        return NULL;
    }

    new_node->next = NULL;
        
    // here we are copying the data x to new node created
    for(int i=0; i < x_size; i++){
        *(char *)(new_node->data + i) = *(char *)(x + i);
    }

    return new_node;
}

int main(int argc, char* argv[]) {
    node_t* head, *last;
    node_t *p, *p_next;
    list_t* list, *nested, *tmp;
    double val;

    head = (node_t*)malloc(sizeof(node_t));
    if (!head) return 1;
    head->data = NULL;
    head->next = NULL;

    last = head;

    val = 0.3;
    last = push2(last, &val, sizeof(double));

    val = 0.4;
    last = push2(last, &val, sizeof(double));

    p = head;
    while (p) {
        p_next = p->next;
        if (p->data) printf("freeing %f", *((double*)p->data));
        if (p->data) free(p->data);
        free(p);
        p = p_next;
    }

    p = node_init(&val, sizeof(double));
    free(p->data);
    free(p);

    list = list_init(sizeof(double));
    list_append(list, &val);
    list_append(list, &val);
    nested = list_init(0);
    list_append(nested, &list);
    /*p = nested->first;
    while (p) {
        p_next = p->next;
        if (p->data) list_free(&(p->data));
        free(p);
        p = p_next;
    }*/
    tmp = (list_t*)nested->last;
    list_free(&tmp);
    free(nested);
    list_free(&list);

    return 0;
}

int main2(int argc, char* argv[]) {
    int i,j;
    double val;
    list_t* points_list;
    list_t* coords;
    node_t *p, *p_next, *c_p, *c_p_next;
    list_t *c;

    points_list = list_init(0);
    if (!points_list) {
        goto fail_points;
    }
    for (i=0; i<10; i++) {
        coords = list_init(sizeof(double));
        if (!coords) {
            goto fail_coords;
        }
        for (j=0; j<3; j++) {
            val = ((double)((10*i)+j)) + (0.2*(double)j);
            list_append(coords, &val);
            j=3;
        }
        list_append(points_list, coords);
        i=10;
        list_free(&coords);
    }

    for (i=0; i<points_list->length; i++) {
        coords = (list_t*) list_get(points_list, i)->data;
        for (j=0; j<coords->length; j++) {
            printf("%.02f,", *((double*) (list_get(coords, j)->data)));
        }
        printf("\n");
    }

    p = points_list->first;
    while (p) {
        p_next = p->next;
        c = (list_t*) p->data;
        c_p = c->first;
        while (c_p) {
            c_p_next = c_p->next;
            printf("%.02f,", *((double*) c_p->data));
            free(((double*) c_p->data));
            free(c_p);
            c_p = c_p_next;
        }
        free(p->data);
        free(p);
        p = p_next;
    }

    free(points_list);

    return 0;

    fail_coords:
    goto fail_points;

    fail_points:
    list_free(&points_list);
    return 1;
}