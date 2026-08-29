/**
 * @file transform.c
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Definitions for spherical distance transforms
 */

#include "transform.h"


#ifdef _WIN32
    #ifndef PYINIT_LIBPLATERECIPY_TRANSFORM
        #define PYINIT_LIBPLATERECIPY_TRANSFORM
        void PyInit_libplaterecipy_transform() {};
    #endif
#endif


// ~~~~~ Callback mechanism for the Python logger ~~~~~

void set_transform_h_logger(
    transform_h_log_func func
) {
    transform_h_logger = func;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~~~~~~~~~~ Internal functions ~~~~~~~~~~~~~~~~~

int add_if_not_stored(int * array, int current_size, int item) {
    bool found = false;
    for (int i = current_size-1; i >= 0; i--) {
        if (array[i] == item) {
            found = true;
            break;
        }
    }
    if (!found) {
        array[current_size++] = item;
    }
    return current_size;
}

bool is_false_border_node(Map * map, int idx, bool * markers, bool safe_mode) {
    if (markers[map->nodes[idx].npy_idx]) {
        // discarding the node if True
        return false;
    }

    // checking for the first-order and second-order neighbours (if safe_mode)
    for (int i = 0; i < map->nodes[idx].num_neighs; i++) {
        const int neigh_idx = map->nodes[idx].neighs[i];

        if (markers[map->nodes[neigh_idx].npy_idx]) {
            // this means the node has a 1st-order True neighbour
            return true;
        }
        
        if (safe_mode) {
            // checking for 2nd neighbours
            for (int j = 0; j < map->nodes[neigh_idx].num_neighs; j++) {
                const int neigh_neigh_idx = map->nodes[neigh_idx].neighs[j]; 

                if (markers[map->nodes[neigh_neigh_idx].npy_idx]) {
                    // this means the node has a 2nd-order True neighbour
                    return true;
                }
            }
        }
        
    }

    return false;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~~~~~~~~~~ Exposed functions ~~~~~~~~~~~~~~~~~

void label_markers_from_map(
    Map *       map,
    bool *      markers,
    int32_t *   labels
) {
    if (transform_h_logger) transform_h_logger("Assigning integer labels to connected patches ...");

    int max_stack_size = 0;
    for (int i = 0; i < map->num_nodes; i++) {
        if (markers[i]) {
            max_stack_size++;
            labels[i] = -1; // marked but not yet processed
        } else {
            labels[i] = 0; // unmarked, so they are assigned zero
        }
    }
    
    int * idx_stack = malloc(max_stack_size * sizeof(int));
    int current_label = 1;

    for (int idx = 0; idx < map->num_nodes; idx++) {
        if (labels[map->nodes[idx].npy_idx] < 0) {
            // an unprocessed marked node
            
            
            // initializing a new stack
            int top = 0;

            // pushing
            idx_stack[top++] = idx;
            labels[map->nodes[idx].npy_idx] = current_label;

            while (top > 0) {
                // popping the stack
                const int popped_idx = idx_stack[--top];

                for (int j = 0; j < map->nodes[popped_idx].num_neighs; j++) {
                    const int neigh_idx = map->nodes[popped_idx].neighs[j];
                    if (labels[map->nodes[neigh_idx].npy_idx] < 0) {
                        // an uprocessed marked neighbour
                        
                        // pushing
                        idx_stack[top++] = neigh_idx;
                        labels[map->nodes[neigh_idx].npy_idx] = current_label;
                    }
                }
            }

            current_label++;
        }
    }

    free(idx_stack);
    if (transform_h_logger) transform_h_logger("... done.");
}


#if defined(__APPLE__) && defined(__MACH__) && !defined(_OPENMP)

struct inverted_fused_distance_threshold_transform_on_map_args {
    Map *       map;
    double *    xs;
    double *    ys;
    double *    zs;
    bool *      arr;
    bool *      arr_out;
    int32_t *   border_nodes_npy_idx;
    int32_t     num_border_nodes;
    double      THRESHOLD_D_SQUARED;
    int32_t     num_threads;
    int32_t     i_thread;
};

void * inverted_fused_distance_threshold_transform_on_map_func (
    void * args
) {
    struct inverted_fused_distance_threshold_transform_on_map_args * targs = 
        (struct inverted_fused_distance_threshold_transform_on_map_args *) args;
    
    int idx_min = (targs->map->num_nodes / targs->num_threads)*(targs->i_thread);
    int idx_max = (targs->map->num_nodes / targs->num_threads)*(targs->i_thread +1);

    // to prevent missing out due to division round downs
    if (targs->i_thread == targs->num_threads - 1) {
        idx_max = targs->map->num_nodes;
    }

    for (int idx = idx_min; idx < idx_max; idx++) {
        const int ref_npy_idx = targs->map->nodes[idx].npy_idx;
        if (targs->arr[ref_npy_idx]) {
            for (int j = 0; j < targs->num_border_nodes; j++) {
                const int border_npy_idx = targs->border_nodes_npy_idx[j];
                
                const double dx = targs->xs[border_npy_idx] - targs->xs[ref_npy_idx];
                const double dy = targs->ys[border_npy_idx] - targs->ys[ref_npy_idx];
                const double dz = targs->zs[border_npy_idx] - targs->zs[ref_npy_idx];
                const double d2 = dx*dx + dy*dy + dz*dz;
                
                if (d2 < targs->THRESHOLD_D_SQUARED) {
                    targs->arr_out[ref_npy_idx] = false;
                    break;
                }
            }
        }
    }
}

void inverted_fused_distance_threshold_transform_on_map(
    Map *       map,
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    double      R,
    double      threshold,
    bool *      arr_out,
    int32_t     num_threads
) {
    if (transform_h_logger) transform_h_logger("Performing an inverted fused distace transfrom ...");
    /**
    THIS IS UNLIKE TYPICAL DISTANCE TRANSFORMS
    THE TRANSFORM IS APPLIED ON THE FALSE VALUES INSTEAD!
    MEANING -> THE FALSE BANDS GET WIDER
    */

    int num_border_nodes = 0;
    int * border_nodes_npy_idx = (int *) malloc(map->num_nodes * sizeof(int));

    
    for (int idx = 0; idx < map->num_nodes; idx++) {
        const int ref_npy_idx = map->nodes[idx].npy_idx;

        // first, initializing output to be the same as input
        arr_out[ref_npy_idx] = arr[ref_npy_idx];
        
        // we keep a list of nodes that comprise the boundary to the FALSE
        // region, so that we only check for proximity to them (instead of all FALSE values)
        if (is_false_border_node(map, idx, arr, true)) {
            border_nodes_npy_idx[num_border_nodes++] = ref_npy_idx;
        }
    }

    

    // if planar, THRESHOLD_D_SQUARED should be merely (R*threshold)^2
    if (map->grid_type == GTYPE_PLANAR) {
        if (transform_h_logger) transform_h_logger("... rescaling R for the distance transform so that R*separation_tolerance == Euclidean distance");
        // we scale the local R so that the constant THRESHOLD_D_SQUARED
        // becomes what we want, instead of the great-circle distance
        R *= threshold*sqrt(1/(2*(1. - cos(threshold))));
    }
    const double THRESHOLD_D_SQUARED = 2.*R*R*(1. - cos(threshold));

    
    if (transform_h_logger) transform_h_logger("... using manual threads to parallelize the loop");
    
    pthread_t tids[num_threads];
    struct inverted_fused_distance_threshold_transform_on_map_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].map                    = map;
        targss[i_thread].xs                     = xs;
        targss[i_thread].ys                     = ys;
        targss[i_thread].zs                     = zs;
        targss[i_thread].arr                    = arr;
        targss[i_thread].arr_out                = arr_out;
        targss[i_thread].border_nodes_npy_idx   = border_nodes_npy_idx;
        targss[i_thread].num_border_nodes       = num_border_nodes;
        targss[i_thread].THRESHOLD_D_SQUARED    = THRESHOLD_D_SQUARED;
        targss[i_thread].num_threads            = num_threads;
        targss[i_thread].i_thread               = i_thread;

        // creating threads
        pthread_create(
            &tids[i_thread], 
            NULL, 
            inverted_fused_distance_threshold_transform_on_map_func, 
            &targss[i_thread]
        );
    }

    // waiting for all threads to be done
    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        pthread_join(tids[i_thread], NULL);
    }
    
    free(border_nodes_npy_idx);
    if (transform_h_logger) transform_h_logger("... done.");
}

#else

void inverted_fused_distance_threshold_transform_on_map(
    Map *       map,
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    double      R,
    double      threshold,
    bool *      arr_out,
    int32_t     num_threads
) {
    if (transform_h_logger) transform_h_logger("Performing an inverted fused distace transfrom ...");
    /**
    THIS IS UNLIKE TYPICAL DISTANCE TRANSFORMS
    THE TRANSFORM IS APPLIED ON THE FALSE VALUES INSTEAD!
    MEANING -> THE FALSE BANDS GET WIDER
    */

    int num_border_nodes = 0;
    int * border_nodes_npy_idx = (int *) malloc(map->num_nodes * sizeof(int));

    
    for (int idx = 0; idx < map->num_nodes; idx++) {
        const int ref_npy_idx = map->nodes[idx].npy_idx;

        // first, initializing output to be the same as input
        arr_out[ref_npy_idx] = arr[ref_npy_idx];
        
        // we keep a list of nodes that comprise the boundary to the FALSE
        // region, so that we only check for proximity to them (instead of all FALSE values)
        if (is_false_border_node(map, idx, arr, true)) {
            border_nodes_npy_idx[num_border_nodes++] = ref_npy_idx;
        }
    }

    

    // if planar, THRESHOLD_D_SQUARED should be merely (R*threshold)^2
    if (map->grid_type == GTYPE_PLANAR) {
        if (transform_h_logger) transform_h_logger("... rescaling R for the distance transform so that R*separation_tolerance == Euclidean distance");
        // we scale the local R so that the constant THRESHOLD_D_SQUARED
        // becomes what we want, instead of the great-circle distance
        R *= threshold*sqrt(1/(2*(1. - cos(threshold))));
    }
    const double THRESHOLD_D_SQUARED = 2.*R*R*(1. - cos(threshold));

    #ifdef _OPENMP
    if (transform_h_logger) transform_h_logger("... using OpenMP to parallelize the loop");
    omp_set_num_threads(num_threads); 

    #pragma omp parallel for shared(map,arr,arr_out,border_nodes_npy_idx,xs,ys,zs)
    #endif
    for (int idx = 0; idx < map->num_nodes; idx++) {
        const int ref_npy_idx = map->nodes[idx].npy_idx;
        if (arr[ref_npy_idx]) {
            for (int j = 0; j < num_border_nodes; j++) {
                const int border_npy_idx = border_nodes_npy_idx[j];
                
                const double dx = xs[border_npy_idx] - xs[ref_npy_idx];
                const double dy = ys[border_npy_idx] - ys[ref_npy_idx];
                const double dz = zs[border_npy_idx] - zs[ref_npy_idx];
                const double d2 = dx*dx + dy*dy + dz*dz;
                
                if (d2 < THRESHOLD_D_SQUARED) {
                    arr_out[ref_npy_idx] = false;
                    break;
                }
            }
        }
    }
    
    free(border_nodes_npy_idx);
    if (transform_h_logger) transform_h_logger("... done.");
}

#endif