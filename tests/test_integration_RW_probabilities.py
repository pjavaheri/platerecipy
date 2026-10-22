"""
@file test_integration_RW_probabilities.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief A script for testing the RW solution probabilities.
"""

import numpy as np
import matplotlib.pyplot as plt
import platerecipy.segmentation as ps
from platerecipy.grid import platerecipy_clib_grid, PlanarGrid, GTYPE_PLANAR
import scipy.ndimage as ndi
from platerecipy import _FLOAT, _INT
import ctypes

## Setup
img = np.array([
    [0.4, 0.2, 0.2],
    [1.0, 0.3, 0.1],
    [0.1, 0.2, 0.4]
])
markers = img < 0.2
labels = markers.copy().astype(_INT)
labels[markers==True] = [1, 2]
beta = 1.

## Part 1 (a basic RW solver)
indices = np.linspace(0, 8, 9, dtype=int).reshape((3, 3))
marked_idx = indices[markers]
unmarked_idx = indices[~markers]
ord2org = np.hstack([marked_idx , unmarked_idx])
org2ord = np.array([np.where(ord2org == i)[0][0] for i in range(ord2org.size)])
L = np.zeros((9, 9))
neighs = [
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
]
for k_ord in range(9):
    k_org = ord2org[k_ord]
    i_org, j_org = k_org//3, k_org%3
    for l in neighs:
        di, dj = l
        ni_org = i_org+di
        nj_org = j_org+dj
        if (0 <= ni_org) and (ni_org < 3) and (0 <= nj_org) and (nj_org  < 3):
            print((i_org, j_org), (ni_org, nj_org))
            delta = img[ni_org, nj_org] - img[i_org,j_org]
            w = np.exp(-delta*delta*beta)
            nk_ord = org2ord[ni_org*3+nj_org]
            L[k_ord, nk_ord] = -w

for d in range(9):
    L[d,d] = -np.sum(L[d, :])

L_U = L[2:, 2:]
B_T = L[2:, :2]

P1_ord = np.hstack([np.array([1., 0.]), np.linalg.solve(L_U, -B_T@np.array([1., 0.]))])
P1_org = np.array(
    [P1_ord[org2ord[i]] for i in range(9)]
).reshape((3,3))
P2_ord = np.hstack([np.array([0., 1.]), np.linalg.solve(L_U, -B_T@np.array([0., 1.]))])
P2_org = np.array(
    [P2_ord[org2ord[i]] for i in range(9)]
).reshape((3,3))

expected_P1 = P1_org
expected_P2 = P2_org
expected_A  = L_U
expected_B  = np.vstack([(-B_T@np.array([1., 0.])).T, (-B_T@np.array([0., 1.])).T]).T


## Part 2 (platerecipy implementation)

map = platerecipy_clib_grid.alloc_rect_map(
    ctypes.c_int32(3),
    ctypes.c_int32(3),
    ctypes.c_int32(GTYPE_PLANAR)
)
# assigning indices 
npy_idxs = np.linspace(0, 8, 9, dtype=_INT)
platerecipy_clib_grid.set_npy_idx(
    map,
    npy_idxs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
)

class TestGrid():
    def __init__(self, map):
        self.map = map
        self._map = map
        self.r = 1.
        self._r = 1.

grid = TestGrid(map)

num_labelled = np.sum(markers)
n_i, n_j = img.shape
A, B, num_labelled = ps.get_AB_from_map(
    data            = img,
    labels          = labels,
    beta            = beta,
    largest_label   = 2,
    grid            = grid
)

X = ps.solve_AX_B(A, B, solver='LU')

X       = X.astype(dtype=_FLOAT, order='C', copy=False)
IDs     = np.zeros((n_i, n_j), dtype=_INT, order='C')
probs   = np.zeros((n_i, n_j, X.shape[1]), _FLOAT, order='C')

ps.platerecipy_clib_segmentation.get_IDs_and_probs_from_X_and_map(
    grid.map,
    X.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ctypes.c_int32(X.shape[0]),
    ctypes.c_int32(X.shape[1]),
    labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
    ctypes.c_int32(num_labelled),
    IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
    probs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
) 


## Analysis
print("\n\n\n\n\n\n\n\n ~~~~~~~~~<< Test Begins >>~~~~~~~~~")

print("\n\nChecking if the LHS is as expected ...")
print(A.todense() - expected_A)
assert np.allclose(A.todense(), expected_A, rtol=1e-06, atol=1e-14), "Expected above ^^^ to be approximately zero. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the LHS is as expected ...")
print(B - expected_B)
assert np.allclose(B, expected_B, rtol=1e-06, atol=1e-14), "Expected above ^^^ to be approximately zero. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the full probability matrix is correctly reproduced ...")
print(probs[:, :, 0]- expected_P1)
assert np.allclose(probs[:, :, 0], expected_P1, rtol=1e-06, atol=1e-14), "Expected above ^^^ to be approximately zero. Integration test failed (x)"
print('\n', probs[:, :, 1]- expected_P2)
assert np.allclose(probs[:, :, 1], expected_P2, rtol=1e-06, atol=1e-14), "Expected above ^^^ to be approximately zero. Integration test failed (x)"
print("    ... integration test passed!")

print("~~~~~~~~~<< Done >>~~~~~~~~~")


fig, axes = plt.subplots(4, 3, figsize=(8, 8))

ax = axes[0][0]
ax.set_title("Input image")
fig.colorbar(ax.imshow(img, vmin=0, vmax=1))

ax = axes[0][1]
ax.set_title("Markers")
fig.colorbar(ax.imshow(markers))

ax = axes[0][2]
ax.set_title("Labels")
fig.colorbar(ax.imshow(labels))

ax = axes[1][0]
ax.set_title("Expected $L_U$")
fig.colorbar(ax.imshow(L_U))

ax = axes[1][1]
ax.set_title("Generated $L_U$")
fig.colorbar(ax.imshow(A.todense()))

ax = axes[1][2]
ax.set_title("Difference")
fig.colorbar(ax.imshow(A.todense()-L_U))

ax = axes[2][0]
ax.set_title("Expected $P_1$")
fig.colorbar(ax.imshow(expected_P1, vmin=0, vmax=1))

ax = axes[2][1]
ax.set_title("Generated $P_1$")
fig.colorbar(ax.imshow(probs[:, :, 0], vmin=0, vmax=1))

ax = axes[2][2]
ax.set_title("Difference")
fig.colorbar(ax.imshow(probs[:, :, 0] - expected_P1))

ax = axes[3][0]
ax.set_title("Expected $P_2$")
fig.colorbar(ax.imshow(expected_P2, vmin=0, vmax=1))

ax = axes[3][1]
ax.set_title("Generated $P_2$")
fig.colorbar(ax.imshow(probs[:, :, 1], vmin=0, vmax=1))

ax = axes[3][2]
ax.set_title("Difference")
fig.colorbar(ax.imshow(probs[:, :, 1] - expected_P2))

fig.tight_layout()
fig.savefig("test_integration_RW_probabilities.png", dpi=600)