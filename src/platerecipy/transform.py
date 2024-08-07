import ctypes
import sysconfig
from pathlib import Path
import numpy as np

# on Linux/Unix, pip should be in the parent directory
_clib_transform = Path(__file__).parent.parent  \
    / Path( "clib_transform" + str(sysconfig.get_config_var('EXT_SUFFIX')) )

"""
A Python access to `clib/transform.h` module.
"""
_clib_transform_mod = ctypes.CDLL(_clib_transform)

def fused_distance_threshold_transform(
    arr: np.ndarray,
    threshold: float,
    double_precision=True
) -> np.ndarray:
    """
    Apply a fused distance threshold transform. The input 2D array `arr` should 
    be a boolean array with cells representing a uniformly spaced grid of a 
    mercator projection of a sphere. The transformation preserves the `True` and
    modifies `False` according to whether the cell is within the vicinity of a
    `True` determined by the great circle angle of separation
    """

    arr_out = np.zeros(arr.shape, type=np.bool)

    if double_precision:
        _clib_transform_mod.fused_distance_threshold_transform_64bit(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int64(arr.shape[0]),
            ctypes.c_int64(arr.shape[1]),
            ctypes.c_double(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
        )
    else:
        _clib_transform_mod.fused_distance_threshold_transform_32bit(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
            ctypes.c_int32(arr.shape[0]),
            ctypes.c_int32(arr.shape[1]),
            ctypes.c_float(threshold),
            arr_out.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
        )