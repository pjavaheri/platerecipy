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
 * Calculates the angle of separation on the great circle separating two points
 * at `(lon1, lat1)` and `(lon2, lat2)` in `float32`.
 * 
 * The equation used is numerically stable for all angles:
 * https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulae
 * 
 * @param lon1 longitude of point 1
 * @param lat1 latitude of point 1
 * @param lon2 longitude of point 2
 * @param lat2 latitude of point 2
 * 
 * @returns full angle of separation in radians
 * 
 * @warning Latitudes are measured from the equator with the domain `[-pi/2, pi/2]`.
 * 
 */
float great_circle_angle_32bit(
    float lon1,
    float lat1,
    float lon2,
    float lat2
);

/**
 * Calculates the angle of separation on the great circle separating two points
 * at `(lon1, lat1)` and `(lon2, lat2)` in `float64`.
 * 
 * The equation used is numerically stable for all angles:
 * https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulae
 * 
 * @param lon1 longitude of point 1
 * @param lat1 latitude of point 1
 * @param lon2 longitude of point 2
 * @param lat2 latitude of point 2
 * 
 * @returns full angle of separation in radians
 * 
 * @warning Latitudes are measured from the equator with the domain `[-pi/2, pi/2]`.
 * 
 */
double great_circle_angle_64bit(
    double lon1,
    double lat1,
    double lon2,
    double lat2
);

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