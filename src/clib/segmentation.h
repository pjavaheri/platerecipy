/**
 * @file transform.h
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Header file for image segmentation
 */

#ifndef CLIB_SEGMENTATION_H
#define CLIB_SEGMENTATION_H

#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>

#ifdef _WIN32
    #define CLIB_EXPORT __declspec(dllexport)
#else
    #define CLIB_EXPORT
#endif

#include "grid.h"

/**
Callback mechanism for communicating with python logger
*/
typedef void (*segmentation_h_log_func)(const char *);
static segmentation_h_log_func segmentation_h_logger = NULL;
CLIB_EXPORT void set_segmentation_h_logger(segmentation_h_log_func func);

// Pi
const double PI     = 3.14159265358979323846;
// Pi / 2
const double PI_2   = 1.57079632679489661923;
// Pi / 4
const double PI_4   = 0.78539816339744830962;
// A big floating point number
const double BIG    = 1e6;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
//                         Random Walker routines                             //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


/**
 * Populates vectors `rows`, `column`, and `values` (to later use it to define 
 * the sparse Laplacian matrix), with values obtained from a Gaussian dependence 
 * on a uniform rectilinear grid. 
 * 
 * @param map a pointer to an initialized Map struct
 * @param image a numpy array of input image
 * @param labels a numpy array of labels
 * @param beta Gaussian scaling factor
 * @param rows_ord vector indicating the ordered row
 * @param columns_ord vector indicating the ordered column
 * @param values vector carrying the corresponding value
 *
 * @returns number of nodes that are labelled
 * 
 * @warning `rows_ord`, `columns_ord`, and `values` must be allocated array of size: 
 *           2*(map->num_edges)
 * @warning arrays `rows_ord`, `columns_ord`, and `values` will be updated.
 */
CLIB_EXPORT int32_t get_ordered_Laplacian_from_map(
    Map *       map,
    double *    image,
    int32_t *   labels,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values
);


/**
 * Populates `M` which is the matrix on the right-hand-side that applies 
 * unit test potential to the segments.
 * 
 * @param map a pointer to an initialized Map struct
 * @param labels a numpy array of labels
 * @param num_labelled number of labelled nodes
 * @param largest_label largest label (also equal to the number of segments)
 * @param M the matrix applying unit test potentials
 */
CLIB_EXPORT void get_ordered_boundary_matrix_from_map(
    Map *       map,
    int32_t *   labels,
    int32_t     largest_label,
    int32_t     num_labelled,
    double *    M
);


/**
 * Obtains plate IDs (`IDs`) and ID probabilities (`probs`) from the ordered 
 * solution of the linear system. 
 * 
 * @param map a pointer to an initialized Map struct
 * @param X 2D array of solutions
 * @param X_i_max number of unmarked nodes
 * @param X_j_max number of segments
 * @param labels a numpy array of labels
 * @param num_labelled number of labelled nodes
 * @param IDs array of IDs
 * @param probs array of probabilities 
 */
CLIB_EXPORT void get_IDs_and_probs_from_X_and_map(
    Map *       map,
    double *    X,
    int32_t     X_i_max,
    int32_t     X_j_max,
    int32_t *   labels,
    int32_t     num_labelled,
    int32_t *   IDs,
    double *    probs
);


#endif