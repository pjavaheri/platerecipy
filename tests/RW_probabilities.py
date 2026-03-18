"""
@file RW_probabilities.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief A script for testing the RW solution probabilities.
"""

import numpy as np
import platerecipy.segmentation as ps
from platerecipy import _FLOAT, _INT
import scipy.ndimage as ndi
import ctypes

# creating the input image
img = np.matrix([
    [0.4, 0.2, 0.3, 0.2],
    [1.0, 0.3, 0.4, 0.1],
    [0.3, 0.8, 0.2, 0.3],
    [0.1, 0.2, 0.9, 0.4]
])

# image shape
n_i, n_j = img.shape

# specifying and integer labelling the markers
markers = img < 0.2
labels, _ = ndi.label(markers)
num_labelled = np.sum(markers)


beta = 2.
A, B, ord2org = ps.get_AB(
    data            = img,
    labels          = labels,
    beta            = beta ,
    num_labelled    = 2,
    largest_label   = 2,
    is_spherical    = False
)

X = ps.solve_AX_B(A, B)

X       = X.astype(dtype=_FLOAT, order='C', copy=False)
IDs     = np.zeros((n_i, n_j), dtype=_INT, order='C')
probs   = np.zeros((n_i, n_j, X.shape[1]), _FLOAT, order='C')

ps.platerecipy_clib_segmentation.get_IDs_and_probs_from_X(
    X.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ctypes.c_int32(X.shape[0]),
    ctypes.c_int32(X.shape[1]),
    labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
    ord2org.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
    ctypes.c_int32(n_i),
    ctypes.c_int32(n_j),
    ctypes.c_int32(num_labelled),
    IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
    probs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
) 


expected_P1 = np.array([
    [0.63477274, 0.69343836, 0.78938582, 0.89469291],
    [0.5235145 , 0.65274011, 0.78002618, 1.        ],
    [0.24999473, 0.46122762, 0.70620423, 0.80514888],
    [0.        , 0.24082412, 0.58393568, 0.72058967]
])

expected_P2 = np.array([
    [0.36522726, 0.30656164, 0.21061418, 0.10530709],
    [0.4764855 , 0.34725989, 0.21997382, 0.        ],
    [0.75000527, 0.53877238, 0.29379577, 0.19485112],
    [1.        , 0.75917588, 0.41606432, 0.27941033]
])

assert np.all(np.abs(probs[:, :, 0] - expected_P1) < 1e-7) and np.all(np.abs(probs[:, :, 1] - expected_P2) < 1e-7)

print('SUCCESS')
