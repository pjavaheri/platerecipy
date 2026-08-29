/**
 * @file segmentation.c
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Definitions for image segmentation
 */

#include "segmentation.h"

#ifdef _WIN32
    #ifndef PYINIT_LIBPLATERECIPY_SEGMENTATION
        #define PYINIT_LIBPLATERECIPY_SEGMENTATION
        void PyInit_libplaterecipy_segmentation() {};
    #endif
#endif

void set_segmentation_h_logger(segmentation_h_log_func func) {
    segmentation_h_logger = func;
}

/**
 * (internal)
 * 
 * Find the index that stores the largest value.
 * 
 * @param a array of values
 * @param i_max size of the array
 * 
 * @returns index
 */
int argmax(
    double * a,
    int i_max
) {
    int index = 0;
    double maxval = a[0];

    for (int i = 1; i < i_max; i++) {
        if (a[i] > maxval) {
            maxval = a[i];
            index = i;
        }
    }

    return index;
}

/**
 * (internal)
 * 
 * Compare two floating point values within a tolerance.
 * 
 * @param a first number
 * @param b second number
 * @param eps acceptable tolerance
 * 
 * @returns index
 */
inline bool approx_equal(double a, double b, double eps) {
    return fabs(a-b) < eps;
}

int32_t get_ordered_Laplacian_from_map(
    Map *       map,
    double *    image,
    int32_t *   labels,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values
) {
    segmentation_h_logger("Entered get_ordered_Laplacian_from_map ...");
    // ~~~~~ assigning ord_idx to nodes ~~~~~
    // first the marked
    int k = 0;
    for (int i = 0; i < map->num_nodes; i++) {
        if (labels[map->nodes[i].npy_idx] > 0) {
            map->nodes[i].ord_idx = k;
            k++;
        }
    }
    const int32_t num_labelled = k;
    // then the unmarked
    for (int i = 0; i < map->num_nodes; i++) {
        if (labels[map->nodes[i].npy_idx] == 0) {
            map->nodes[i].ord_idx = k;
            k++;    
        }
    }
    segmentation_h_logger("... done with ord_idx assignment");
    // ~~~~~ constructing the ordered Laplacian ~~~~~
    k = 0;
    for (int i = 0; i < map->num_nodes; i++) {
        for (int j = 0; j < map->nodes[i].num_neighs; j++) {
            const Node * node = map->nodes + (i);
            const Node * neigh = map->nodes + (map->nodes[i].neighs[j]);
            rows_ord[k] = node->ord_idx;
            columns_ord[k] = neigh->ord_idx;
            const double contrast = image[node->npy_idx] - image[neigh->npy_idx];
            double metric_correction = 1./map->nodes[i].edge_lengths[j];
            if (((isinf(metric_correction)) || (isnan(metric_correction))) 
                || (metric_correction > BIG)) {
                metric_correction = BIG;
            }
            values[k] = -exp(-beta*(contrast*contrast))*metric_correction;
            // the negative comes from -wij in the definition of L_ij for
            // adjacent nodes
            k++;
        }
    }

    segmentation_h_logger("... done with get_ordered_Laplacian_from_map.");
    return num_labelled;
}


void get_ordered_boundary_matrix_from_map(
    Map *       map,
    int32_t *   labels,
    int32_t     largest_label,
    int32_t     num_labelled,
    double *    M
) {
    segmentation_h_logger("Entered get_ordered_boundary_matrix_from_map ...");
    for (int i = 0; i < map->num_nodes; i++) {
        const int label_i = labels[map->nodes[i].npy_idx];
        const int ord_idx = map->nodes[i].ord_idx;
        if (label_i > 0) {
            if (ord_idx >= num_labelled) {
                segmentation_h_logger("Something bad has happned 294");
            }         
            for (int s = 0; s < largest_label; s++) {
                if (label_i == s+1) {
                    M[ord_idx*largest_label + s] = 1.;
                }
            }
        }
    }
    segmentation_h_logger("... done with get_ordered_boundary_matrix_from_map.");
}


void get_IDs_and_probs_from_X_and_map(
    Map *       map,
    double *    X,
    int32_t     X_i_max,
    int32_t     X_j_max,
    int32_t *   labels,
    int32_t     num_labelled,
    int32_t *   IDs,
    double *    probs
) {
    for (int i = 0; i < map->num_nodes; i++) {
        const int npy_idx = map->nodes[i].npy_idx;
        if (labels[npy_idx] > 0) {
            IDs[npy_idx] = labels[npy_idx];
            probs[npy_idx*X_j_max + labels[npy_idx]-1] = 1.;
        } else {
            IDs[npy_idx] = argmax(
                X + (map->nodes[i].ord_idx-num_labelled)*X_j_max,
                X_j_max
            ) + 1;
            for (int j = 0; j < X_j_max; j++) {
                probs[npy_idx*X_j_max + j] = X[(map->nodes[i].ord_idx-num_labelled)*X_j_max + j];
            }
        }
    }
}


/**
 * Under development
 */
int get_segment_boundaries(
    int32_t *   IDs,
    int32_t     n_i,
    int32_t     n_j,
    bool        is_wraparound,
    bool *      boundaries
){
    int i;
    int j;
    for (i = 1; i < n_i-1; i++) {
        for (j = 1; j < n_j-1; j++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                          + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }
    }
    i = 1;
    for (j = 1; j < n_j-1; j++) {
        //const int ref_ID = IDs;
        const int sum = IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                      + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);
    }
    i = n_i-1;
    for (j = 1; j < n_j-1; j++) {
        //const int ref_ID = IDs;
        const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                      + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);
    }
    if (is_wraparound) {
        int sum;
        i = 0;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        i = 0;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j-1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        i = n_i-1;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        i = n_i-1;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j-1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        j = 0;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (n_j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                          + IDs[(i  )*n_j + (n_j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (n_j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }

        j = n_j-1;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (0  )]
                          + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (0  )]
                          + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (0  )];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }
    } else {
        int sum;
        i = 0;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)]
            + IDs[(i  )*n_j + (n_j-1)] + IDs[(i+1)*n_j + (n_j-1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        i = 0;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j-1)]
            + IDs[(i  )*n_j + (0  )] + IDs[(i+1)*n_j + (0  )];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        i = n_i-1;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
            + IDs[(i  )*n_j + (n_j-1)] + IDs[(i-1)*n_j + (n_j-1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        i = n_i-1;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j-1)]
            + IDs[(i  )*n_j + (0  )] + IDs[(i-1)*n_j + (0  )];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        j = 0;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                          + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }

        j = n_j-1;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] 
                          + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] 
                          + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }
        
    }
    return 0;
}