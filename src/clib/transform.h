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
 * @warning `arr_out` must be initialized with `-1` for the plate of interest 
 *          and `-2` for all other plates/regions.
 */
void single_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
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
 * @param arr_out for the output distances
 */
void full_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t *   plate_IDs,
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
 * @warning `arr_out` must be initialized with `-1` for the plate of interest 
 *          and `-2` for all other plates/regions.
 */
void single_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t     num_points,
    double      R,
    double *    arr_out,
    int64_t     num_threads
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
 */
void full_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t *   plate_IDs,
    int64_t     num_points,
    double      R,
    double *    arr_out,
    int64_t     num_threads
);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                    LEGACY CODE
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * Spherical distance and threshold transforms on a mercator projection. This 
 * function updates `arr_out` such that all `true` in `arr` is preserved and
 * the `false` depending on whether their great circle angular distance from a 
 * local `true` is less than `threshold` becomes a `true`.
 * 
 * @param arr a pointer to a 2D array of shape `i_max` by `j_max`
 * @param i_max extent of `arr` along the first index
 * @param j_max extent of `arr` along the second index
 * @param threshold angle of separation threshold in radians
 * @param arr_out a pointer to a 2D array for an output
 * 
 * @warning This function is intended to be called from a Python interface.
 * @warning Both `arr` and `arr_out` are considered to be row-major arrays.
 * @warning `i` and `j` in `arr` and `arr_out` are uniformly spaced longitudes
 *          and latitudes. It is assumed that extent is `[0,2pi]x[0,pi]`.
 * @warning `arr_out` must be initialized with `false`.
 */
void sph_fused_distance_threshold_transform_32bit(
    bool * arr,
    int32_t i_max,
    int32_t j_max,
    float threshold,
    bool * arr_out
);

/**
 * Spherical distance and threshold transforms on a mercator projection. This 
 * function updates `arr_out` such that all `true` in `arr` is preserved and
 * the `false` depending on whether their great circle angular distance from a 
 * local `true` is less than `threshold` becomes a `true`. This function is
 * threaded to help with the processing time.
 * 
 * @param arr a pointer to a 2D array of shape `i_max` by `j_max`
 * @param i_max extent of `arr` along the first index
 * @param j_max extent of `arr` along the second index
 * @param threshold angle of separation threshold in radians
 * @param arr_out a pointer to a 2D array for an output
 * @param num_threads number of threads
 * 
 * @warning This function is intended to be called from a Python interface.
 * @warning Both `arr` and `arr_out` are considered to be row-major arrays.
 * @warning `i` and `j` in `arr` and `arr_out` are uniformly spaced longitudes
 *          and latitudes. It is assumed that extent is `[0,2pi]x[0,pi]`.
 * @warning `arr_out` must be initialized with `false`.
 */
void sph_fused_distance_threshold_transform_32bit_threaded(
    bool * arr,
    int32_t i_max,
    int32_t j_max,
    float threshold,
    bool * arr_out,
    int32_t num_threads
);

/**
 * Spherical distance and threshold transforms on a mercator projection. This 
 * function updates `arr_out` such that all `true` in `arr` is preserved and
 * the `false` depending on whether their great circle angular distance from a 
 * local `true` is less than `threshold` becomes a `true`. This version is 
 * double precision.
 * 
 * @param arr a pointer to a 2D array of shape `i_max` by `j_max`
 * @param i_max extent of `arr` along the first index
 * @param j_max extent of `arr` along the second index
 * @param threshold angle of separation threshold in radians
 * @param arr_out a pointer to a 2D array for an output
 * 
 * @warning This function is intended to be called from a Python interface.
 * @warning Both `arr` and `arr_out` are considered to be row-major arrays.
 * @warning `i` and `j` in `arr` and `arr_out` are uniformly spaced longitudes
 *          and latitudes. It is assumed that extent is `[0,2pi]x[0,pi]`.
 * @warning `arr_out` must be initialized with `false`.
 */
void sph_fused_distance_threshold_transform_64bit(
    bool * arr,
    int64_t i_max,
    int64_t j_max,
    double threshold,
    bool * arr_out
);

/**
 * Spherical distance and threshold transforms on a mercator projection. This 
 * function updates `arr_out` such that all `true` in `arr` is preserved and
 * the `false` depending on whether their great circle angular distance from a 
 * local `true` is less than `threshold` becomes a `true`. This version is 
 * double precision. This function is threaded to help with the processing time.
 * 
 * @param arr a pointer to a 2D array of shape `i_max` by `j_max`
 * @param i_max extent of `arr` along the first index
 * @param j_max extent of `arr` along the second index
 * @param threshold angle of separation threshold in radians
 * @param arr_out a pointer to a 2D array for an output
 * @param num_threads number of threads
 * 
 * @warning This function is intended to be called from a Python interface.
 * @warning Both `arr` and `arr_out` are considered to be row-major arrays.
 * @warning `i` and `j` in `arr` and `arr_out` are uniformly spaced longitudes
 *          and latitudes. It is assumed that extent is `[0,2pi]x[0,pi]`.
 * @warning `arr_out` must be initialized with `false`.
 */
void sph_fused_distance_threshold_transform_64bit_threaded(
    bool * arr,
    int64_t i_max,
    int64_t j_max,
    double threshold,
    bool * arr_out,
    int64_t num_threads
);

#endif