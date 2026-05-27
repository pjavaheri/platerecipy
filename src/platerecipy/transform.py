"""
@file transform.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Internal module for transformation functions.
"""

import logging
log = logging.getLogger(__name__)

import os
import ctypes
import numpy as np

from . import _INT, _FLOAT, _IS_WINDOWS

# importing shared libraries
import sysconfig

module_path         = os.path.dirname(os.path.abspath(__file__))
module_path         = os.path.abspath(os.path.join(module_path, os.pardir))
shared_object_path  = os.path.join(
    module_path, 
    "libplaterecipy_transform" + sysconfig.get_config_var('EXT_SUFFIX')
)

log.debug("Loading libplaterecipy_transform shared library")
"""
A Python access to `clib/transform.h` module.
"""
platerecipy_clib_transform = ctypes.CDLL(shared_object_path)

def single_plate_interior_distance_transform(
    xs                  : np.ndarray,
    ys                  : np.ndarray,
    zs                  : np.ndarray,
    plate_indicators    : np.ndarray,
    R                   = None,
    num_threads         = 1
) -> np.ndarray:
    """
    Performs a spherical (i.e., great-circle) distance transform on a given 
    patch/mask indicated by `plate_indicators`.

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
    
    num_threads : int, optional
        Number of threads to divide the task. It is possible instead pass the
        string `'auto'` which results in using number of threads equal to the
        number of CPUs.
        default = 1
    
    Returns
    -------
        np.ndarray

    """
    if _IS_WINDOWS:
        log.info("Multithreaded functionalities are not available on Windows and `num_threads` will be ignored.")
        num_threads = 1

    if R is None:
        R = np.mean(np.sqrt(np.pow(xs**2 + ys**2 + zs**2, 0.5)))

    if num_threads == 'auto':
        num_threads = os.cpu_count()

    # to interface with C
    xs = xs.astype(dtype=_FLOAT, order='C', copy=False)
    ys = ys.astype(dtype=_FLOAT, order='C', copy=False)
    zs = zs.astype(dtype=_FLOAT, order='C', copy=False)
    plate_indicators = plate_indicators.astype(order='C', copy=False)
    
    arr_out = -1. * np.ones(plate_indicators.shape, dtype=_FLOAT, order='C')
    arr_out[~plate_indicators] = -2.

    if num_threads > 1:
        platerecipy_clib_transform.single_plate_interior_distance_transform_64bit_threaded(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(arr_out.size),
            ctypes.c_double(R),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(num_threads)
        )
    else:
        platerecipy_clib_transform.single_plate_interior_distance_transform_64bit(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(arr_out.size),
            ctypes.c_double(R),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        )
    
    return arr_out


def full_plate_interior_distance_transform(
    xs          : np.ndarray,
    ys          : np.ndarray,
    zs          : np.ndarray,
    plate_IDs   : np.ndarray,
    R           = None,
    num_threads = 1
) -> np.ndarray:
    """
    Performs a spherical (i.e., great-circle) distance transform on every patch
    marked by integer labels in `plate_IDs`.

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
    
    num_threads : int, optional
        Number of threads to divide the task. It is possible instead pass the
        string `'auto'` which results in using number of threads equal to the
        number of CPUs.
        default = 1
    
    Returns
    -------
        np.ndarray
    """
    if _IS_WINDOWS:
        log.info("Multithreaded functionalities are not available on Windows and `num_threads` will be ignored.")
        num_threads = 1

    if R is None:
        R = np.mean(np.sqrt(xs**2 + ys**2 + zs**2))

    if num_threads == 'auto':
        num_threads = os.cpu_count()

    # to interface with C
    xs = xs.astype(dtype=_FLOAT, order='C', copy=False)
    ys = ys.astype(dtype=_FLOAT, order='C', copy=False)
    zs = zs.astype(dtype=_FLOAT, order='C', copy=False)
    plate_IDs = plate_IDs.astype(dtype=_INT, order='C', copy=False)

    arr_out = -1. * np.ones(plate_IDs.shape, dtype=_FLOAT, order='C')

    if num_threads > 1:
        platerecipy_clib_transform.full_plate_interior_distance_transform_64bit_threaded(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            plate_IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            ctypes.c_int32(arr_out.size),
            ctypes.c_double(R),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(num_threads)
        )
    else:
        platerecipy_clib_transform.full_plate_interior_distance_transform_64bit(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            plate_IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            ctypes.c_int32(arr_out.size),
            ctypes.c_double(R),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        )

    return arr_out


def fused_distance_threshold_transform(
    xs          : np.ndarray,
    ys          : np.ndarray,
    zs          : np.ndarray,
    arr         : np.ndarray,
    threshold   : float,
    R           = None
) -> np.ndarray:
    """
    Apply a fused distance threshold transform. The transformation preserves the 
    `True` and modifies `False` according to whether the cell is within the 
    vicinity of a `True` determined by the great circle angle of separation.

    Parameters
    ----------
    xs : np.ndarray
        Cartesian x positions.
    
    ys : np.ndarray
        Cartesian y positions.
    
    zs : np.ndarray
        Cartesian z positions.

    arr : np.ndarray
        Input boolean array to apply the transformation.
    
    threshold : float
        Angular threshold (great circle separation) in radians.
    
    R : float, optional
        Reference radius for calculating arc angles. If not provided, average of 
        `(xs**2 + ys**2 + zs**2)**0.5` will be used.
    
    
    Returns
    -------
        np.ndarray

    """
    if R is None:
        R = np.mean(np.sqrt(xs**2 + ys**2 + zs**2))

    # to interface with C
    xs = xs.astype(dtype=_FLOAT, order='C', copy=False)
    ys = ys.astype(dtype=_FLOAT, order='C', copy=False)
    zs = zs.astype(dtype=_FLOAT, order='C', copy=False)
    arr = arr.astype(dtype=bool, order='C', copy=False)
    
    arr_out = np.zeros(arr.shape, dtype=bool, order='C')

    platerecipy_clib_transform.fused_distance_threshold_transform_64bit(
        xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
        ctypes.c_int32(arr.shape[0]),
        ctypes.c_int32(arr.shape[1]),
        ctypes.c_double(R),
        ctypes.c_double(threshold),
        arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    )
        
    return arr_out


def gridded_fused_distance_threshold_transform(
    xs          : np.ndarray,
    ys          : np.ndarray,
    zs          : np.ndarray,
    arr         : np.ndarray,
    threshold   : float,
    R           = None,
    num_threads = 1
) -> np.ndarray:
    """
    Apply a fused distance threshold transform. The input 2D array `arr` should 
    be a boolean array with cells representing a uniformly spaced grid of a 
    mercator projection of a sphere. The transformation preserves the `True` and
    modifies `False` according to whether the cell is within the vicinity of a
    `True` determined by the great circle angle of separation

    Parameters
    ----------
    xs : np.ndarray
        Cartesian x positions.
    
    ys : np.ndarray
        Cartesian y positions.
    
    zs : np.ndarray
        Cartesian z positions.

    arr : np.ndarray
        Input boolean array to apply the transformation.
    
    threshold : float
        Angular threshold (great circle separation) in radians.
    
    R : float, optional
        Reference radius for calculating arc angles. If not provided, average of 
        `(xs**2 + ys**2 + zs**2)**0.5` will be used.
    
    num_threads : int, optional
        Number of threads to divide the task. It is possible instead pass the
        string `'auto'` which results in using number of threads equal to the
        number of CPUs.
        default = 1
    
    
    Returns
    -------
        np.ndarray

    """
    if _IS_WINDOWS:
        log.info("Multithreaded functionalities are not available on Windows and `num_threads` will be ignored.")
        num_threads = 1

    if R is None:
        log.debug("No radius provided. Estimating it from the node coordinates.")
        R = np.mean(np.sqrt(xs**2 + ys**2 + zs**2))

    if num_threads == 'auto':
        num_threads = os.cpu_count()

    # to interface with C
    xs = xs.astype(dtype=_FLOAT, order='C', copy=False)
    ys = ys.astype(dtype=_FLOAT, order='C', copy=False)
    zs = zs.astype(dtype=_FLOAT, order='C', copy=False)
    arr = arr.astype(dtype=bool, order='C', copy=False)
    
    arr_out = np.zeros(arr.shape, dtype=bool, order='C')

    if num_threads > 1:
        log.debug("Calling gridded_fused_distance_threshold_transform_64bit_threaded from C on %d threads", num_threads)
        platerecipy_clib_transform.gridded_fused_distance_threshold_transform_64bit_threaded(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(arr.shape[0]),
            ctypes.c_int32(arr.shape[1]),
            ctypes.c_double(R),
            ctypes.c_double(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(num_threads),
        )
    else:
        log.debug("Calling gridded_fused_distance_threshold_transform_64bit from C")
        platerecipy_clib_transform.gridded_fused_distance_threshold_transform_64bit(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(arr.shape[0]),
            ctypes.c_int32(arr.shape[1]),
            ctypes.c_double(R),
            ctypes.c_double(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
        )
        
    return arr_out


def pgridded_fused_distance_threshold_transform(
    xs          : np.ndarray,
    ys          : np.ndarray,
    zs          : np.ndarray,
    theta_range : tuple,
    phi_range   : tuple,
    arr         : np.ndarray,
    threshold   : float,
    R           = None,
    num_threads = 1
) -> np.ndarray:
    """
    Apply a fused distance threshold transform. The input 2D array `arr` should 
    be a boolean array with cells representing a uniformly spaced grid of a 
    mercator projection of a sphere. The transformation preserves the `True` and
    modifies `False` according to whether the cell is within the vicinity of a
    `True` determined by the great circle angle of separation. This function 
    operates on a partial sphere defined by `theta_range` and `phi_range`.

    Parameters
    ----------
    xs : np.ndarray
        Cartesian x positions.
    
    ys : np.ndarray
        Cartesian y positions.
    
    zs : np.ndarray
        Cartesian z positions.
    
    theta_range : tuple,
        Range of the polar angle.
    
    phi_range : tuple,
        Range of the azimuthal angle.

    arr : np.ndarray
        Input boolean array to apply the transformation.
    
    threshold : float
        Angular threshold (great circle separation) in radians.
    
    R : float, optional
        Reference radius for calculating arc angles. If not provided, average of 
        `(xs**2 + ys**2 + zs**2)**0.5` will be used.
    
    num_threads : int, optional
        Number of threads to divide the task. It is possible instead pass the
        string `'auto'` which results in using number of threads equal to the
        number of CPUs.
        default = 1
    
    
    Returns
    -------
        np.ndarray

    """
    if _IS_WINDOWS:
        log.info("Multithreaded functionalities are not available on Windows and `num_threads` will be ignored.")
        num_threads = 1
        
    if R is None:
        log.debug("No radius provided. Estimating it from the node coordinates.")
        R = np.mean(np.sqrt(xs**2 + ys**2 + zs**2))

    if num_threads == 'auto':
        num_threads = os.cpu_count()

    # to interface with C
    xs = xs.astype(dtype=_FLOAT, order='C', copy=False)
    ys = ys.astype(dtype=_FLOAT, order='C', copy=False)
    zs = zs.astype(dtype=_FLOAT, order='C', copy=False)
    arr = arr.astype(dtype=bool, order='C', copy=False)
    
    arr_out = np.zeros(arr.shape, dtype=bool, order='C')

    if num_threads > 1:
        log.debug("Calling pgridded_fused_distance_threshold_transform_64bit_threaded from C on %d threads", num_threads)
        platerecipy_clib_transform.pgridded_fused_distance_threshold_transform_64bit_threaded(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(arr.shape[0]),
            ctypes.c_int32(arr.shape[1]),
            ctypes.c_double(R),
            ctypes.c_double(theta_range[0]),
            ctypes.c_double(theta_range[1]),
            ctypes.c_double(phi_range[0]),
            ctypes.c_double(phi_range[1]),
            ctypes.c_double(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(num_threads),
        )
    else:
        log.debug("Calling pgridded_fused_distance_threshold_transform_64bit from C")
        platerecipy_clib_transform.pgridded_fused_distance_threshold_transform_64bit(
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(arr.shape[0]),
            ctypes.c_int32(arr.shape[1]),
            ctypes.c_double(R),
            ctypes.c_double(theta_range[0]),
            ctypes.c_double(theta_range[1]),
            ctypes.c_double(phi_range[0]),
            ctypes.c_double(phi_range[1]),
            ctypes.c_double(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
        )
        
    return arr_out
