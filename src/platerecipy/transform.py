"""
@file transform.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for transformation functions.
"""
import os
import ctypes
import sysconfig
from pathlib import Path
import numpy as np

# on Linux/Unix, pip should be in the parent directory
_platerecipy_clib_transform = Path(__file__).parent.parent  \
    / Path( "platerecipy_clib_transform" + str(sysconfig.get_config_var('EXT_SUFFIX')) )

"""
A Python access to `clib/transform.h` module.
"""
_platerecipy_clib_transform_mod = ctypes.CDLL(_platerecipy_clib_transform)

def sph_distance_transform(
    arr: np.ndarray,
    threshold: float,
    double_precision=True,
    num_threads=1
) -> np.ndarray:
    """
    Apply a distance transform. The input 2D array `arr` should be a boolean 
    array with cells representing a uniformly spaced grid of a mercator 
    projection of a sphere. The transformation preserves the `True` as 0 
    otherwise `False` cells will be recorded as the distance to the nearest cell 
    determined by the great circle angle of separation.

    Parameters
    ----------
    arr : np.ndarray
        Input boolean array to apply the transformation.
    
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

    if double_precision:
        arr_out = -1. * np.ones(arr.shape, dtype=np.float64)
        if num_threads > 1:
            _platerecipy_clib_transform_mod.sph_distance_transform_64bit_threaded(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(arr.shape[0]),
                ctypes.c_int64(arr.shape[1]),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int64(num_threads)
            )
        else:
            _platerecipy_clib_transform_mod.sph_distance_transform_64bit(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(arr.shape[0]),
                ctypes.c_int64(arr.shape[1]),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            )
    else:
        arr_out = -1. * np.ones(arr.shape, dtype=np.float32)
        if num_threads > 1:
            _platerecipy_clib_transform_mod.sph_distance_transform_32bit_threaded(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(arr.shape[0]),
                ctypes.c_int32(arr.shape[1]),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_int32(num_threads)
            )
        else:
            _platerecipy_clib_transform_mod.sph_distance_transform_32bit(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(arr.shape[0]),
                ctypes.c_int32(arr.shape[1]),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )

    return arr_out

def sph_fused_distance_threshold_transform(
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
            _platerecipy_clib_transform_mod.sph_fused_distance_threshold_transform_64bit_threaded(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(arr.shape[0]),
                ctypes.c_int64(arr.shape[1]),
                ctypes.c_double(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(num_threads)
            )
        else:
            _platerecipy_clib_transform_mod.sph_fused_distance_threshold_transform_64bit(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int64(arr.shape[0]),
                ctypes.c_int64(arr.shape[1]),
                ctypes.c_double(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
            )
    else:
        if num_threads > 1:
            _platerecipy_clib_transform_mod.sph_fused_distance_threshold_transform_32bit_threaded(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(arr.shape[0]),
                ctypes.c_int32(arr.shape[1]),
                ctypes.c_float(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(num_threads)
            )
        else:
            _platerecipy_clib_transform_mod.sph_fused_distance_threshold_transform_32bit(
                arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                ctypes.c_int32(arr.shape[0]),
                ctypes.c_int32(arr.shape[1]),
                ctypes.c_float(threshold),
                arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
            )

    return arr_out