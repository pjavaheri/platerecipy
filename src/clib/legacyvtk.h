/**
 * @file legacyvtk.h
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Header file for legacy VTK output
 */

#ifndef CLIB_LEGACYVTK_H
#define CLIB_LEGACYVTK_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#if defined(_WIN32)
    #define CLIB_EXPORT __declspec(dllexport)
#else
    #define CLIB_EXPORT
#endif

/**
 * Initializes an ASCII VTK file and sets up a mesh grid (points + faces).
 * 
 * @note This function can be used for both a rectangular grid and a partial
 * spherical grid.
 * 
 * @param adr string address of the vtk output file
 * @param xs array of Cartesian x coordinates
 * @param ys array of Cartesian y coordinates
 * @param zs array of Cartesian z coordinates
 * @param i_max extent along the first dimension
 * @param j_max extent along the second dimension
 */
CLIB_EXPORT int make_rectangular_vtk_grid(
    char *      adr,
    double *    xs,
    double *    ys,
    double *    zs,
    int32_t     i_max,
    int32_t     j_max
);

/**
 * Adds a float field data to an already initialized ASCII VTK file.
 * 
 * @note This function can be used for both a rectangular grid and a partial
 * spherical grid.
 * 
 * @param adr string address of the vtk output file
 * @param field_name string name of the field
 * @param field array of field data
 * @param i_max extent along the first dimension
 * @param j_max extent along the second dimension
 */
CLIB_EXPORT int add_rectangular_vtk_float_field(
    char *      adr,
    char *      field_name,
    double *    field,
    int32_t     i_max,
    int32_t     j_max
);

/**
 * Adds an int field data to an already initialized ASCII VTK file.
 * 
 * @note This function can be used for both a rectangular grid and a partial
 * spherical grid.
 * 
 * @param adr string address of the vtk output file
 * @param field_name string name of the field
 * @param field array of field data
 * @param i_max extent along the first dimension
 * @param j_max extent along the second dimension
 */
CLIB_EXPORT int add_rectangular_vtk_int_field(
    char *      adr,
    char *      field_name,
    int32_t *   field,
    int32_t     i_max,
    int32_t     j_max
);

/**
 * Initializes an ASCII VTK file and sets up a mesh grid (points + faces).
 * 
 * @note This function can only be used for a spherical grid.
 * @warning The input data must be structured so that the first and last row 
 * correspond to the north and south poles, and the last column is equal to the 
 * first column (thus ignored).
 * 
 * @param adr string address of the vtk output file
 * @param xs array of Cartesian x coordinates
 * @param ys array of Cartesian y coordinates
 * @param zs array of Cartesian z coordinates
 * @param i_max extent along the first dimension
 * @param j_max extent along the second dimension
 */
CLIB_EXPORT int make_spherical_vtk_grid(
    char *      adr,
    double *    xs,
    double *    ys,
    double *    zs,
    int32_t     i_max,
    int32_t     j_max
);

/**
 * Adds field data to an already initialized ASCII VTK file.
 * 
 * @note This function can only be used for a spherical grid.
 * @warning The input data must be structured so that the first and last row 
 * correspond to the north and south poles, and the last column is equal to the 
 * first column (thus ignored).
 * 
 * @param adr string address of the vtk output file
 * @param field_name string name of the field
 * @param field array of field data
 * @param i_max extent along the first dimension
 * @param j_max extent along the second dimension
 */
CLIB_EXPORT int add_spherical_vtk_float_field(
    char *      adr,
    char *      field_name,
    double *    field,
    int32_t     i_max,
    int32_t     j_max
);

/**
 * Adds field data to an already initialized ASCII VTK file.
 * 
 * @note This function can only be used for a spherical grid.
 * @warning The input data must be structured so that the first and last row 
 * correspond to the north and south poles, and the last column is equal to the 
 * first column (thus ignored).
 * 
 * @param adr string address of the vtk output file
 * @param field_name string name of the field
 * @param field array of field data
 * @param i_max extent along the first dimension
 * @param j_max extent along the second dimension
 */
CLIB_EXPORT int add_spherical_vtk_int_field(
    char *      adr,
    char *      field_name,
    int32_t *   field,
    int32_t     i_max,
    int32_t     j_max
);

#endif