/**
 * @file transform.h
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Header file for image transformations
 */

#ifndef CLIB_TRANSFORM_H
#define CLIB_TRANSFORM_H

#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>

// Pi
const double PI     = 3.14159265358979323846;
// Pi / 2
const double PI_2   = 1.57079632679489661923;
// Pi / 4
const double PI_4   = 0.78539816339744830962;

// Boolean true
const int TRUE      = 1;
// Boolean false
const int FALSE     = 0;

/**
 * Calculate the angle of separation on the great circle separating two points
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

#endif