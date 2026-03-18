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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>


// Pi
const double PI     = 3.14159265358979323846;
// Pi / 2
const double PI_2   = 1.57079632679489661923;
// Pi / 4
const double PI_4   = 0.78539816339744830962;
// A big floating point number
const double BIG    = 1e6;

/** NOTE
 * Spherical version of each function (if applicable) has a trailing `_sph` in 
 * its name. 
 */

/** NOTE
 * Functions under the `Interface` banner are intended to be called from outside.
 */

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
//                         Random Walker routines                             //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


/**
 * Populates the `edges` array with the vertices denoting corresponding 
 * edges (connections) for a rectilinear grid.
 * 
 * @param edges array pointer
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * 
 * @warning `edges` must be an allocated array of size: 
 *           2*( (n_i-1)*n_j + n_i*(n_j-1) ) 
 * @warning `edges` will be updated
 */
int populate_edges(
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j
);

/**
 * Populates the `edges` array with the vertices denoting corresponding 
 * edges (connections) for a spherical grid.
 * 
 * The fundamental difference between this function and its rectilinear 
 * counter part (`populate_edges`) is that the last element of the first
 * row (i.e., [0, n_j-1]) and the first element of the last row (i.e., 
 * [n_i-1, 0]) are considered as the north and south poles with their 
 * respective connections. 
 * 
 * @param edges array pointer
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * 
 * @warning `edges` must be an allocated array of size: 
 *           2*( (2*n_i-3)*n_j ) 
 * @warning `edges` will be updated
 */
int populate_edges_sph(
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j
);


/**
 * Populates vectors `rows`, `column`, and `values` (to later use it to define 
 * the sparse Laplacian matrix), with values obtained from a Gaussian dependence 
 * on a uniform rectilinear grid. 
 * 
 * @param image a 2D array of input image
 * @param edges array of edges
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param beta Gaussian scaling factor
 * @param rows vector indicating the row
 * @param columns vector indicating the column
 * @param values vector carrying the corresponding value
 * 
 * @warning `rows`, `column`, and `values` must be allocated array of size: 
 *           2*( (n_i-1)*n_j + n_i*(n_j-1) ) 
 * @warning arrays `rows`, `columns`, and `values` will be updated.
 */
int get_Laplacian_from_edges(
    double *    image,
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    int32_t *   rows,
    int32_t *   columns,
    double *    values
);


/**
 * Populates vectors `rows`, `column`, and `values` (to later use it to define 
 * the sparse Laplacian matrix), with values obtained from a Gaussian dependence 
 * on a slice of a uniform spherical grid with an additional multiplicative 
 * metric correction factor to make up for the polar distortion. The metric 
 * correction is capped at constant `BIG`.
 * 
 * @param image a 2D array of input image
 * @param edges array of edges
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param theta_min lower range of theta
 * @param theta_max upper range of theta
 * @param phi_min lower range of phi
 * @param phi_max upper range of phi
 * @param beta Gaussian scaling factor
 * @param rows vector indicating the row
 * @param columns vector indicating the column
 * @param values vector carrying the corresponding value
 * 
 * @warning `rows`, `column`, and `values` must be allocated array of size: 
 *           2*( (n_i-1)*n_j + n_i*(n_j-1) ) 
 * @warning arrays `rows`, `columns`, and `values` will be updated.
 */
int get_Laplacian_from_edges_psph(
    double *    image,
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j,
    double      theta_min,
    double      theta_max,
    double      phi_min,
    double      phi_max,
    double      beta,
    int32_t *   rows,
    int32_t *   columns,
    double *    values
);


/**
 * Populates vectors `rows`, `column`, and `values` (to later use it to define 
 * the sparse Laplacian matrix), with values obtained from a Gaussian dependence 
 * on a uniform spherical grid with an additional multiplicative metric correction
 * factor to make up for the polar distortion. The metric correction is capped 
 * at constant `BIG`.
 * 
 * @param image a 2D array of input image
 * @param edges array of edges
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param beta Gaussian scaling factor
 * @param rows vector indicating the row
 * @param columns vector indicating the column
 * @param values vector carrying the corresponding value
 * 
 * @warning `rows`, `column`, and `values` must be allocated array of size: 
 *           2*( (2*n_i-3)*n_j ) 
 * @warning arrays `rows`, `columns`, and `values` will be updated.
 */
int get_Laplacian_from_edges_sph(
    double *    image,
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    int32_t *   rows,
    int32_t *   columns,
    double *    values
);


/**
 * Populate `org2ord` and `ord2org` to allow for a conversion between the original 
 * (i.e., row-major ordering from the original image) to the ordered (i.e., sorted
 * by wether they are labelled or not) for a rectilinear grid.
 * 
 * @param labels 2D array of labels 
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param org2ord vector mapping original ordering to sorted ordering
 * @param ord2org vector mapping sorted ordering to original ordering
 * 
 * @warning `ord2org` and `org2ord` must be allocated arrays of size: n_i*n_j
 
 * @warning arrays `org2ord` and `ord2org` will be updated.
 */
int get_original_to_ordered_mapping(
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   org2ord,
    int32_t *   ord2org
);


/**
 * Populate `org2ord` and `ord2org` to allow for a conversion between the original 
 * (i.e., row-major ordering from the original image) to the ordered (i.e., sorted
 * by wether they are labelled or not) for a spherical grid.
 * 
 * @param labels 2D array of labels 
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param org2ord vector mapping original ordering to sorted ordering
 * @param ord2org vector mapping sorted ordering to original ordering
 * 
 * @warning `ord2org` and `org2ord` must be allocated arrays of size:
 *           (n_i-2)*n_j + 2
 * @warning arrays `org2ord` and `ord2org` will be updated.
 * @warning the conversion vectors are shifted on the original ordering as the 
 *          first item is at n_j-1 (i.e., the north pole).
 */
int get_original_to_ordered_mapping_sph(
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   org2ord,
    int32_t *   ord2org
);


/**
 * Populate `rows_ord` and `columns_ord` by taking the original ordering and 
 * sorting them for a rectilinear grid.
 * 
 * @param org2ord vector mapping original ordering to sorted ordering
 * @param rows vector indicating the row
 * @param columns vector indicating the column
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param rows_ord ordered vector indicating the row
 * @param columns_ord ordered vector indicating the column
 * 
 * @warning `rows_ord` and `columns_ord` must be allocated arrays of size: n_i*n_j
 * @warning  arrays `rows_ord` and `columns_ord` will be updated.
 */
int order_Laplacian(
    int32_t *   org2ord,
    int32_t *   rows,
    int32_t *   columns,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   rows_ord,
    int32_t *   columns_ord
);


/**
 * Populate `rows_ord` and `columns_ord` by taking the original ordering and 
 * sorting them for a spherical grid.
 * 
 * @param org2ord vector mapping original ordering to sorted ordering
 * @param rows vector indicating the row
 * @param columns vector indicating the column
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param rows_ord ordered vector indicating the row
 * @param columns_ord ordered vector indicating the column
 * 
 * @warning `rows_ord` and `columns_ord` must be allocated arrays of size:
 *          (n_i-2)*n_j + 2
 * @warning  arrays `rows_ord` and `columns_ord` will be updated.
 * @warning the conversion vectors are to be shifted on the original ordering as the 
 *          first item is at n_j-1 (i.e., the north pole).
 */
int order_Laplacian_sph(
    int32_t *   org2ord,
    int32_t *   rows,
    int32_t *   columns,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   rows_ord,
    int32_t *   columns_ord
);


// ~~~~~~~~~~~~~~~~~~~~~~~~~ Interface functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


/**
 * Populates `rows_ord`, `column_ord`, and `values` to be used as vectors to
 * construct a sparse matrix, with ordered indices.
 * 
 * @warning This function provides an adequate initialization and wrapping of 
 * functions: 
 * `populate_edges()`,
 * `get_Laplacian_from_edges()`,
 * `get_original_to_ordered_mapping()`, and
 * `order_Laplacian()`.
 * 
 * @param data a 2D array of input image
 * @param labels the 2D array of markers
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param beta Gaussian scaling factor
 * @param rows_ord ordered vector indicating the row
 * @param columns_ord ordered vector indicating the column
 * @param values vector carrying the corresponding value
 * @param ord2org vector mapping original ordering to sorted ordering
 */
int get_ordered_Laplacian_vectors(
    double *    data,
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values,
    int32_t *   ord2org
);


/**
 * Populates `rows_ord`, `column_ord`, and `values` to be used as vectors to
 * construct a sparse matrix, with ordered indices.
 * 
 * @warning This function provides an adequate initialization and wrapping of 
 * functions: 
 * `populate_edges()`,
 * `get_Laplacian_from_edges_psph()`,
 * `get_original_to_ordered_mapping()`, and
 * `order_Laplacian()`.
 * 
 * @param data a 2D array of input image
 * @param labels the 2D array of markers
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param theta_min lower range of theta
 * @param theta_max upper range of theta
 * @param phi_min lower range of phi
 * @param phi_max upper range of phi
 * @param beta Gaussian scaling factor
 * @param rows_ord ordered vector indicating the row
 * @param columns_ord ordered vector indicating the column
 * @param values vector carrying the corresponding value
 * @param ord2org vector mapping original ordering to sorted ordering
 */
int get_ordered_Laplacian_vectors_psph(
    double *    data,
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    double      theta_min,
    double      theta_max,
    double      phi_min,
    double      phi_max,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values,
    int32_t *   ord2org
);


/**
 * Populates `rows_ord`, `column_ord`, and `values` to be used as vectors to
 * construct a sparse matrix, with ordered indices, implemented for a spherical
 * input.
 * 
 * @warning Input arrays `data` and `labels` must be trimmed off of their 
 * redundant last column, and shape variables `n_i` and `n_j` must reflect the
 * dimensions of the trimmed array.
 * 
 * @note This function provides an adequate initialization and wrapping of 
 * functions: 
 * `populate_edges()`,
 * `get_Laplacian_from_edges()`,
 * `get_original_to_ordered_mapping()`, and
 * `order_Laplacian()`.
 * 
 * @param data a 2D array of input image
 * @param labels the 2D array of markers
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param beta Gaussian scaling factor
 * @param rows_ord ordered vector indicating the row
 * @param columns_ord ordered vector indicating the column
 * @param values vector carrying the corresponding value
 * @param ord2org vector mapping original ordering to sorted ordering
 */
int get_ordered_Laplacian_vectors_sph(
    double *    data,
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values,
    int32_t *   ord2org
);



/**
 * Populates `M` which is the matrix on the right-hand-side that applies 
 * unit test potential to the segments.
 * 
 * 
 * @param labels the 2D array of markers
 * @param ord2org vector mapping original ordering to sorted ordering
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param num_labelled number of labelled nodes
 * @param largest_label largest label (also equal to the number of segments)
 * @param M the matrix applying unit test potentials
 */
int get_ordered_boundary_matrix(
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t     largest_label,
    double *    M
);


/**
 * Populates `M` which is the matrix on the right-hand-side that applies 
 * unit test potential to the segments.
 * 
 * @warning Input arrays `data` and `labels` must be trimmed off of their 
 * redundant last column, and shape variables `n_i` and `n_j` must reflect the
 * dimensions of the trimmed array.
 * 
 * 
 * @param labels the 2D array of markers
 * @param ord2org vector mapping original ordering to sorted ordering
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction
 * @param num_labelled number of labelled nodes
 * @param largest_label largest label (also equal to the number of segments)
 * @param M the matrix applying unit test potentials
 */
int get_ordered_boundary_matrix_sph(
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t     largest_label,
    double *    M
);


/**
 * Obtains plate IDs (`IDs`) and ID probabilities (`probs`) from the ordered 
 * solution of the linear system. 
 * 
 * @param X 2D array of solutions
 * @param X_i_max number of unmarked nodes
 * @param X_j_max number of segments
 * @param labels 2D array of markers
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction 
 * @param num_labelled number of labelled nodes
 * @param IDs array of IDs
 * @param probs array of probabilities 
 */
int get_IDs_and_probs_from_X(
    double *    X,
    int32_t     X_i_max,
    int32_t     X_j_max,
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t *   IDs,
    double *    probs
);


/**
 * Obtains plate IDs (`IDs`) and ID probabilities (`probs`) from the ordered 
 * solution of the linear system. 
 * 
 * @warning Input arrays `data` and `labels` must be trimmed off of their 
 * redundant last column, and shape variables `n_i` and `n_j` must reflect the
 * dimensions of the trimmed array.
 * 
 * @param X 2D array of solutions
 * @param X_i_max number of unmarked nodes
 * @param X_j_max number of segments
 * @param labels 2D array of markers
 * @param n_i number of grid points along the polar direction
 * @param n_j number of grid points along the azimuthal direction 
 * @param num_labelled number of labelled nodes
 * @param IDs array of IDs
 * @param probs array of probabilities 
 */
int get_IDs_and_probs_from_X_sph(
    double *    X,
    int32_t     X_i_max,
    int32_t     X_j_max,
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t *   IDs,
    double *    probs
);

#endif