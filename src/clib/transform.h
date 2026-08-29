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

#include "grid.h"

#ifdef _WIN32
    #define CLIB_EXPORT __declspec(dllexport)
#else
    #define CLIB_EXPORT
#endif

#ifdef _OPENMP
    #include <omp.h>
#else
    #if defined(__APPLE__) && defined(__MACH__)
        #include <pthread.h>
    #endif
#endif


// ~~~~~ Callback mechanism for the Python logger ~~~~~

typedef void (*transform_h_log_func)(const char *);
static transform_h_log_func transform_h_logger = NULL;
CLIB_EXPORT void set_transform_h_logger(transform_h_log_func func);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



// Pi
const double PI     = 3.14159265358979323846;
// Pi / 2
const double PI_2   = 1.57079632679489661923;
// Pi / 4
const double PI_4   = 0.78539816339744830962;
// A big floating point number
const double BIG    = 1e6;



/**
 * Labels True connected patches by positive integer values 
 * (starting from 1) and leaves False intervening bands as False.
 * 
 * @param map a pointer to a Map struct
 * @param markers a pointer to the boolean input array
 * @param labels a pointer to an allocated array to store the output
 *
 * @warning labels will be overwritten.
 */
CLIB_EXPORT void label_markers_from_map(
    Map *       map,
    bool *      markers,
    int32_t *   labels
);


/**
 * Performs an inverted fused distance transform. This means,
 * all True regions that are close enough (defined by threshold)
 * to False regions will be turned to False in the output.
 *
 * @note R corresponds to the radius of the sphere, so that
 * the threshold will be the great-circle distance in radians.
 * Otherwise, R will effectively non-dimensionalize the distances.
 *
 * @param map a pointer to a Map struct
 * @param xs Cartesian x coordinates
 * @param ys Cartesian y coordinates
 * @param zs Cartesian z coordinates 
 * @param arr input boolean array
 * @param R the radius of the sphere
 * @param threshold distance threshold in radians
 * @param arr_out for the output transformed boolean array
 */
CLIB_EXPORT void inverted_fused_distance_threshold_transform_on_map(
    Map *       map,
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    double      R,
    double      threshold,
    bool *      arr_out,
    int32_t     num_threads
);

#endif