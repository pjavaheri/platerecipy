"""
@file transform.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Internal module for transformation functions.
"""

import logging
log = logging.getLogger(__name__)

import os
import numpy as np

from . import _INT, _FLOAT, _IS_WINDOWS
from .grid import Grid, Map, SphericalGrid

# ~~~~~~~~~~~~~~ importing shared libraries ~~~~~~~~~~~~~
#import weakref
import os
import sysconfig
import ctypes

module_path         = os.path.dirname(os.path.abspath(__file__))
module_path         = os.path.abspath(os.path.join(module_path, os.pardir))
shared_object_path  = os.path.join(
    module_path, 
    "libplaterecipy_transform" + sysconfig.get_config_var('EXT_SUFFIX')
)

log.debug("Loading platerecipy_clib_transform shared library")
"""
A Python access to `clib/transform.h` module.
"""
platerecipy_clib_transform = ctypes.CDLL(shared_object_path)


# platerecipy_clib_transform.label_markers_from_map() function signature
platerecipy_clib_transform.label_markers_from_map.argtypes = [
    ctypes.POINTER(Map),
    ctypes.POINTER(ctypes.c_bool), 
    ctypes.POINTER(ctypes.c_int32)
]
platerecipy_clib_transform.label_markers_from_map.restype = None


# platerecipy_clib_transform.inverted_fused_distance_threshold_transform_on_map() function signature
platerecipy_clib_transform.inverted_fused_distance_threshold_transform_on_map.argtypes = [
    ctypes.POINTER(Map),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool), 
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_bool),
    ctypes.c_int32
]
platerecipy_clib_transform.inverted_fused_distance_threshold_transform_on_map.restype = None


# establishing a callback mechanism for enabling logging from c
callback_func_type = ctypes.CFUNCTYPE(None, ctypes.c_char_p)

def _transform_h_c_log(msg):
    log.debug(f"[platerecipy_clib_transform callback]: {msg.decode()}")

_transform_h_callback = callback_func_type(_transform_h_c_log) # not to be removed
platerecipy_clib_transform.set_transform_h_logger(_transform_h_callback)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def label_markers(
    grid    : Grid,
    markers : np.ndarray,
) -> np.ndarray:
    labels = np.zeros(markers.shape, dtype=_INT)

    log.debug("Calling label_markers_from_map form C")
    platerecipy_clib_transform.label_markers_from_map(
        grid.map,
        markers.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
    )

    if isinstance(grid, SphericalGrid):
        grid.enforce_data_consistency(labels)

    return labels


def distance_transform(
    grid        : Grid,
    markers     : np.ndarray,
    threshold   : float,
    num_threads = 1
)-> np.ndarray:
    out = np.zeros_like(markers)

    log.debug("Calling inverted_fused_distance_threshold_transform_on_map form C")
    platerecipy_clib_transform.inverted_fused_distance_threshold_transform_on_map(
        grid.map,
        grid.xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        grid.ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        grid.zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        markers.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
        ctypes.c_double(grid._r),
        ctypes.c_double(threshold),
        out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
        ctypes.c_int32(num_threads)
    )

    if isinstance(grid, SphericalGrid):
        grid.enforce_data_consistency(out)

    return out

