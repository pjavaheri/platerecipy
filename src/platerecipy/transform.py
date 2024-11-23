"""
@file transform.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for transformation functions.
"""
import os
import ctypes
from pathlib import Path
import numpy as np

# importing shared libraries
import sysconfig

module_path         = os.path.dirname(os.path.abspath(__file__))
module_path         = os.path.abspath(os.path.join(module_path, os.pardir))
shared_object_path  = os.path.join(
    module_path, 
    "platerecipy_clib_transform" + sysconfig.get_config_var('EXT_SUFFIX')
)
print(f"shared_object_path={shared_object_path}")
"""
A Python access to `clib/transform.h` module.
"""
platerecipy_clib_transform = ctypes.CDLL(shared_object_path)

def fused_distance_threshold_transform(
    xs                  : np.ndarray,
    ys                  : np.ndarray,
    zs                  : np.ndarray,
    arr                 : np.ndarray,
    R                   : float,
    threshold           : float,
    double_precision    = True
) -> np.ndarray:
    """
    Apply a fused distance threshold transform. The input 2D array `arr` should 
    be a boolean array with cells representing a uniformly spaced grid of a 
    mercator projection of a sphere. The transformation preserves the `True` and
    modifies `False` according to whether the cell is within the vicinity of a
    `True` determined by the great circle angle of separation

    Parameters
    ----------
    arr : np.ndarray
        Input boolean array to apply the transformation.
    
    threshold : float
        Angular threshold (great circle separation) in radians.
    
    double_precision : bool, default=True
        If the floating point operations must be carried out in double precision.
    
    
    Returns
    -------
        np.ndarray

    """
    print("in fused")
    

    arr_out = np.zeros(arr.shape, dtype=bool)

    if double_precision:
        platerecipy_clib_transform.gridded_fused_distance_threshold_transform_64bit(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int64(arr.shape[0]),
            ctypes.c_int64(arr.shape[1]),
            ctypes.c_double(R),
            ctypes.c_double(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
        )
    else:
        platerecipy_clib_transform.gridded_fused_distance_threshold_transform_64bit(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(arr.shape[0]),
            ctypes.c_int32(arr.shape[1]),
            ctypes.c_float(R),
            ctypes.c_float(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
        )
    print("out fused")
    return arr_out


def legacy_sph_fused_distance_threshold_transform(
    arr: np.ndarray,
    threshold: float,
    double_precision=True,
    num_threads=1
) -> np.ndarray:
    """
    Apply a fused distance threshold transform. The input 2D array `arr` should 
    be a boolean array with cells representing a uniformly spaced grid of a 
    mercator projection of a sphere. The transformation preserves the `True` and
    modifies `False` according to whether the cell is within the vicinity of a
    `True` determined by the great circle angle of separation

    Parameters
    ----------
    arr : np.ndarray
        Input boolean array to apply the transformation.
    
    threshold : float
        Angular threshold (great circle separation) in radians.
    
    double_precision : bool, default=True
        If the floating point operations must be carried out in double precision.
    
    num_threads : int, default=1
        Number of threads to divide the task. It is possible instead pass the
        string `'auto'` which results in using number of threads equal to the
        number of CPUs.
    
    Returns
    -------
        np.ndarray

    Warning
    -------
    Aside from the usual computational trade off between parallel threads and 
    communication, in this function, the algorithm benefits from progression 
    through the array. When the array is divided between threads, the benefit 
    of adjacent threads can get diminished. Although in most cases, using several
    threads (4, 8, or 16) is effective in speeding up calculation. 
    """

    if num_threads == 'auto':
        num_threads = os.cpu_count()

    arr_out = np.zeros(arr.shape, dtype=bool)

    if double_precision:
        if num_threads > 1:
            platerecipy_clib_transform.sph_fused_distance_threshold_transform_64bit_threaded(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(arr.shape[0]),
                ctypes.c_int64(arr.shape[1]),
                ctypes.c_double(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(num_threads)
            )
        else:
            platerecipy_clib_transform.sph_fused_distance_threshold_transform_64bit(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(arr.shape[0]),
                ctypes.c_int64(arr.shape[1]),
                ctypes.c_double(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
            )
    else:
        if num_threads > 1:
            platerecipy_clib_transform.sph_fused_distance_threshold_transform_32bit_threaded(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(arr.shape[0]),
                ctypes.c_int32(arr.shape[1]),
                ctypes.c_float(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(num_threads)
            )
        else:
            platerecipy_clib_transform.sph_fused_distance_threshold_transform_32bit(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(arr.shape[0]),
                ctypes.c_int32(arr.shape[1]),
                ctypes.c_float(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
            )

    return arr_out

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def single_plate_interior_distance_transform(
    xs                  : np.ndarray,
    ys                  : np.ndarray,
    zs                  : np.ndarray,
    plate_indicators    : np.ndarray,
    R                   = None,
    double_precision    = True,
    num_threads         = 1
) -> np.ndarray:
    """
    TBC

    Parameters
    ----------
    xs : np.ndarray
        Cartesian x positions.
    
    ys : np.ndarray
        Cartesian y positions.
    
    zs : np.ndarray
        Cartesian z positions.
    
    plate_indicator : np.ndarray
        Input boolean array to apply the transformation.
    
    R : float, optional
        Reference radius for calculating arc angles. If not provided, average of 
        `(xs**2 + ys**2 + zs**2)**0.5` will be used.
    
    double_precision : bool, optional
        If the floating point operations must be carried out in double precision.
        default = True
    
    num_threads : int, optional
        Number of threads to divide the task. It is possible instead pass the
        string `'auto'` which results in using number of threads equal to the
        number of CPUs.
        default = 1
    
    Returns
    -------
        np.ndarray

    """
    if R is None:
        R = np.mean(np.sqrt(np.pow(xs**2 + ys**2 + zs**2, 0.5)))

    if num_threads == 'auto':
        num_threads = os.cpu_count()

    if double_precision:
        arr_out = -1. * np.ones(plate_indicators.shape, dtype=np.float64)
        arr_out[~plate_indicators] = -2.

        if num_threads > 1:
            platerecipy_clib_transform.single_plate_interior_distance_transform_64bit_threaded(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int64(arr_out.size),
                ctypes.c_double(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int64(num_threads)
            )
        else:
            platerecipy_clib_transform.single_plate_interior_distance_transform_64bit(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int64(arr_out.size),
                ctypes.c_double(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            )
    else:
        arr_out = -1. * np.ones(plate_indicators.shape, dtype=np.float32)
        arr_out[~plate_indicators] = -2.
        if num_threads > 1:
            platerecipy_clib_transform.single_plate_interior_distance_transform_32bit_threaded(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_int32(arr_out.size),
                ctypes.c_float(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_int32(num_threads)
            )
        else:
            platerecipy_clib_transform.single_plate_interior_distance_transform_32bit(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_int32(arr_out.size),
                ctypes.c_float(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )

    return arr_out


def full_plate_interior_distance_transform(
    xs                  : np.ndarray,
    ys                  : np.ndarray,
    zs                  : np.ndarray,
    plate_IDs           : np.ndarray,
    R                   = None,
    double_precision    = True,
    num_threads         = 1
) -> np.ndarray:
    """
    TBC

    Parameters
    ----------
    xs : np.ndarray
        Cartesian x positions.
    
    ys : np.ndarray
        Cartesian y positions.
    
    zs : np.ndarray
        Cartesian z positions.
    
    plate_IDs : np.ndarray
        Plate ID numbers.
    
    R : float, optional
        Reference radius for calculating arc angles. If not provided, average of 
        `(xs**2 + ys**2 + zs**2)**0.5` will be used.
    
    double_precision : bool, optional
        If the floating point operations must be carried out in double precision.
        default = True
    
    num_threads : int, optional
        Number of threads to divide the task. It is possible instead pass the
        string `'auto'` which results in using number of threads equal to the
        number of CPUs.
        default = 1
    
    Returns
    -------
        np.ndarray

    """
    if R is None:
        R = np.mean(np.sqrt(np.pow(xs**2 + ys**2 + zs**2, 0.5)))

    if num_threads == 'auto':
        num_threads = os.cpu_count()

    if double_precision:
        arr_out = -1. * np.ones(plate_IDs.shape, dtype=np.float64)

        if num_threads > 1:
            platerecipy_clib_transform.full_plate_interior_distance_transform_64bit_threaded(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                plate_IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
                ctypes.c_int64(arr_out.size),
                ctypes.c_double(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int64(num_threads)
            )
        else:
            platerecipy_clib_transform.full_plate_interior_distance_transform_64bit(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                plate_IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
                ctypes.c_int64(arr_out.size),
                ctypes.c_double(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            )
    else:
        arr_out = -1. * np.ones(plate_IDs.shape, dtype=np.float32)

        if num_threads > 1:
            platerecipy_clib_transform.full_plate_interior_distance_transform_32bit_threaded(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                plate_IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                ctypes.c_int32(arr_out.size),
                ctypes.c_float(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_int32(num_threads)
            )
        else:
            platerecipy_clib_transform.full_plate_interior_distance_transform_32bit(
                xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                plate_IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                ctypes.c_int32(arr_out.size),
                ctypes.c_float(R),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )

    return arr_out
