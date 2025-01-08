/**
 * @file transform.h
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Header file for image transformations
 */

#ifndef CLIB_TRANSFORM_H
#define CLIB_TRANSFORM_H

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


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//      single_plate_interior_distance_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the spherical distance 
 * transform is to be applied on the interior of only a single plates, indicated 
 * by the initial -1 values for plate interior and -2 values for point not on 
 * the plate.
 */


/**
 * Finds the spherical (geodesic) distance transform for a single plate. 
 * The plate interior is denoted by `arr_out` that is initialized by `-1` (the 
 * plate of interest) and `-2` (all other regions).
 * 
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param num_points length of the array (total number of points)
 * @param R the radius of the sphere
 * @param arr_out for the output distances
 * 
 * @returns 0 if no error
 * 
 * @warning `arr_out` must be initialized with `-1` for the plate of interest 
 *          and `-2` for all other plates/regions.
 */
int single_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t     num_points,
    double      R,
    double *    arr_out
);

/**
 * Finds the spherical (geodesic) distance transform for a single plate. 
 * The plate interior is denoted by `arr_out` that is initialized by `-1` (the 
 * plate of interest) and `-2` (all other regions).
 * 
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param num_points length of the array (total number of points)
 * @param R the radius of the sphere
 * @param arr_out for the output distances
 * @param num_threads number of threads
 * 
 * @returns 0 if no error
 * 
 * @warning `arr_out` must be initialized with `-1` for the plate of interest 
 *          and `-2` for all other plates/regions.
 */
int single_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t     num_points,
    double      R,
    double *    arr_out,
    int64_t     num_threads
);




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       full_plate_interior_distance_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the spherical distance 
 * transform is to be applied on the interior of all plates, indicated in the
 * `plate_IDs` array.
 */


/**
 * Finds the spherical (geodesic) distance transform for all plates. 
 * 
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param plate_IDs plate IDs
 * @param num_points length of the array (total number of points)
 * @param R the radius of the sphere
 * @param arr_out for the output distances
 * 
 * @returns 0 if no error
 */
int full_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t *   plate_IDs,
    int64_t     num_points,
    double      R,
    double *    arr_out
);

/**
 * Finds the spherical (geodesic) distance transform for all plates. 
 * 
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param plate_IDs plate IDs
 * @param num_points length of the array (total number of points)
 * @param R the radius of the sphere
 * @param arr_out for the output distances,
 * @param num_threads number of threads
 * 
 * @returns 0 if no error
 */
int full_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t *   plate_IDs,
    int64_t     num_points,
    double      R,
    double *    arr_out,
    int64_t     num_threads
);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     fused_distance_threshold_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the fused distance 
 * transform is applied on any array of booleans. In all, the original 
 * array `arr` must be different from the output `arr_out`.
 */


/**
 * Performs a distance transform and applies a threshold transform such that 
 * the resulting array contains `true` for all points withing the distance 
 * of `threshold` measured from the great-circle distances on the sphere.
 * 
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param arr initial boolean array
 * @param i_max length of the array along the first dimension
 * @param j_max length of the array along the second dimension
 * @param R the radius of the sphere
 * @param threshold distance threshold in radians
 * @param arr_out for the output transformed boolean array
 * 
 * @returns 0 if no error
 */
int fused_distance_threshold_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int64_t     num_points,
    double      R,
    double      threshold,
    bool *      arr_out
);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     gridded_fused_distance_threshold_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the fused distance 
 * transform is applied on a 2D gridded array of booleans. In all, the original 
 * array `arr` must be different from the output `arr_out`. In all these functions,
 * all arrays are presumed to be gridded as a uniform spherical grid and have 
 * identical shapes.
 */


/**
 * Performs a distance transform and applies a threshold transform such that 
 * the resulting array contains `true` for all points withing the distance 
 * of `threshold` measured from the great-circle distances on the sphere.
 * 
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param arr initial 2D gridded boolean array
 * @param i_max length of the array along the first dimension
 * @param j_max length of the array along the second dimension
 * @param R the radius of the sphere
 * @param threshold distance threshold in radians
 * @param arr_out for the output transformed boolean array
 * 
 * @returns 0 if no error
 */
int gridded_fused_distance_threshold_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int64_t     i_max,
    int64_t     j_max,
    double      R,
    double      threshold,
    bool *      arr_out
);

/**
 * Performs a distance transform and applies a threshold transform such that 
 * the resulting array contains `true` for all points withing the distance 
 * of `threshold` measured from the great-circle distances on the sphere.
 * This is the threaded version of `gridded_fused_distance_threshold_transform_64bit`.
 * 
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param arr initial 2D gridded boolean array
 * @param i_max length of the array along the first dimension
 * @param j_max length of the array along the second dimension
 * @param R the radius of the sphere
 * @param threshold distance threshold in radians
 * @param arr_out for the output transformed boolean array
 * @param num_threads number of threads
 * 
 * @returns 0 if no error
 */
int gridded_fused_distance_threshold_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int64_t     i_max,
    int64_t     j_max,
    double      R,
    double      threshold,
    bool *      arr_out,
    int64_t     num_threads
);

#endif