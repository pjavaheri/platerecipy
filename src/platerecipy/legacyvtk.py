"""
@file legacyvtk.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for legacy VTK output functions.
"""

import logging
log = logging.getLogger(__name__)

import ctypes
import numpy as np

from . import _INT, _FLOAT

# importing shared libraries
import os
import sysconfig

module_path         = os.path.dirname(os.path.abspath(__file__))
module_path         = os.path.abspath(os.path.join(module_path, os.pardir))
shared_object_path  = os.path.join(
    module_path, 
    "platerecipy_clib_legacyvtk" + sysconfig.get_config_var('EXT_SUFFIX')
)

log.debug("Loading platerecipy_clib_legacyvtk.so")
"""
A Python access to `clib/legacyvtk.h` module.
"""
platerecipy_clib_legacyvtk = ctypes.CDLL(shared_object_path)


def make_rectangular_vtk(
    fname: str, 
    xs: np.ndarray, ys: np.ndarray, zs: np.ndarray,
    fields: list,
    field_names: list
):
    N, M = xs.shape
    # to interface with C
    xs = xs.astype(dtype=_FLOAT, order='C', copy=False)
    ys = ys.astype(dtype=_FLOAT, order='C', copy=False)
    zs = zs.astype(dtype=_FLOAT, order='C', copy=False)

    log.debug("Calling make_rectangular_vtk_grid from C")
    platerecipy_clib_legacyvtk.make_rectangular_vtk_grid(
        ctypes.c_char_p(fname.encode('utf-8')),
        xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int32(N),
        ctypes.c_int32(M),
    )

    for k in range(len(fields)):
        if fields[k].dtype == _INT:
            field = fields[k].astype(dtype=_INT, order='C', copy=False)
            log.debug("Calling add_rectangular_vtk_int_field from C")
            platerecipy_clib_legacyvtk.add_rectangular_vtk_int_field(
                ctypes.c_char_p(fname.encode('utf-8')),
                ctypes.c_char_p(field_names[k].encode('utf-8')),
                field.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                ctypes.c_int32(N),
                ctypes.c_int32(M),
            )
        else:
            field = fields[k].astype(dtype=_FLOAT, order='C', copy=False)
            log.debug("Calling add_rectangular_vtk_float_field from C")
            platerecipy_clib_legacyvtk.add_rectangular_vtk_float_field(
                ctypes.c_char_p(fname.encode('utf-8')),
                ctypes.c_char_p(field_names[k].encode('utf-8')),
                field.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int32(N),
                ctypes.c_int32(M),
            )

        
    

def make_spherical_vtk(
    fname: str, 
    xs: np.ndarray, ys: np.ndarray, zs: np.ndarray,
    fields: list,
    field_names: list
):
    N, M = xs.shape
    # to interface with C
    xs = xs.astype(dtype=_FLOAT, order='C', copy=False)
    ys = ys.astype(dtype=_FLOAT, order='C', copy=False)
    zs = zs.astype(dtype=_FLOAT, order='C', copy=False)

    log.debug("Calling make_spherical_vtk_grid from C")
    platerecipy_clib_legacyvtk.make_spherical_vtk_grid(
        ctypes.c_char_p(fname.encode('utf-8')),
        xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int32(N),
        ctypes.c_int32(M),
    )

    for k in range(len(fields)):
        if fields[k].dtype == _INT:
            field = fields[k].astype(dtype=_INT, order='C', copy=False)
            log.debug("Calling add_spherical_vtk_int_field from C")
            platerecipy_clib_legacyvtk.add_spherical_vtk_int_field(
                ctypes.c_char_p(fname.encode('utf-8')),
                ctypes.c_char_p(field_names[k].encode('utf-8')),
                field.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                ctypes.c_int32(N),
                ctypes.c_int32(M),
            )
        else:
            field = fields[k].astype(dtype=_FLOAT, order='C', copy=False)
            log.debug("Calling add_spherical_vtk_float_field from C")
            platerecipy_clib_legacyvtk.add_spherical_vtk_float_field(
                ctypes.c_char_p(fname.encode('utf-8')),
                ctypes.c_char_p(field_names[k].encode('utf-8')),
                field.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int32(N),
                ctypes.c_int32(M),
            )
    
    