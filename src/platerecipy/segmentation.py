"""
@file segmentation.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Internal module for segmentation functions.
"""
import logging
log = logging.getLogger(__name__)

import ctypes
import numpy as np
import scipy as sp

from . import _INT, _FLOAT
from .grid import Grid, PartialSphericalGrid, SphericalGrid, Map

# ~~~~~~~~~~~~~~ importing shared libraries ~~~~~~~~~~~~~
import os
import sysconfig

module_path         = os.path.dirname(os.path.abspath(__file__))
module_path         = os.path.abspath(os.path.join(module_path, os.pardir))
shared_object_path  = os.path.join(
    module_path, 
    "libplaterecipy_segmentation" + sysconfig.get_config_var('EXT_SUFFIX')
)

log.debug("Loading libplaterecipy_segmentation shared library")
"""
A Python access to `clib/segmentation.h` module.
"""
platerecipy_clib_segmentation = ctypes.CDLL(shared_object_path)


# platerecipy_clib_segmentation.get_ordered_Laplacian_from_map() function signature
platerecipy_clib_segmentation.get_ordered_Laplacian_from_map.argtypes = [
    ctypes.POINTER(Map),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.POINTER(ctypes.c_int32),
    ctypes.POINTER(ctypes.c_double)
]
platerecipy_clib_segmentation.get_ordered_Laplacian_from_map.restype = ctypes.c_int32


# platerecipy_clib_segmentation.get_ordered_boundary_matrix_from_map() function signature
platerecipy_clib_segmentation.get_ordered_boundary_matrix_from_map.argtypes = [
    ctypes.POINTER(Map),
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_double)
]
platerecipy_clib_segmentation.get_ordered_boundary_matrix_from_map.restype = None


# platerecipy_clib_segmentation.get_IDs_and_probs_from_X_and_map() function signature
platerecipy_clib_segmentation.get_IDs_and_probs_from_X_and_map.argtypes = [
    ctypes.POINTER(Map),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.POINTER(ctypes.c_double)
]
platerecipy_clib_segmentation.get_IDs_and_probs_from_X_and_map.restype = None


# establishing a callback mechanism for enabling logging from c
callback_func_type = ctypes.CFUNCTYPE(None, ctypes.c_char_p)

def _segmentation_h_c_log(msg):
    log.debug(f"[platerecipy_clib_segmentation callback]: {msg.decode()}")

_segmentation_h_callback = callback_func_type(_segmentation_h_c_log) # not to be removed
platerecipy_clib_segmentation.set_segmentation_h_logger(_segmentation_h_callback)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_AB_from_map(
    data             : np.ndarray,
    labels           : np.ndarray,
    beta             : float,
    largest_label    : int,
    grid             : Grid
) -> tuple:
    """
    Constructs the RW linear system by generating the right-hand-side, `'A'`, and
    left-hand-side, `'B'`, matrices.

    Parameters
    ----------
    data : np.ndarray,
        Input 2D grid.

    labels : np.ndarray,
        Labels (markers) corresponding to `data`.

    beta : float,
        The beta factor for connection weights.

    largest_label : int,
        The largest label (i.e., the number of segments).

    grid : Grid,
        An instance of the grid object to perform segmentation on.
    
    Returns
    -------
    tuple(`A`, `RHS`)


    Warning
    -------
    If `grid_type == sph`, it is assumed that the first axis 
    (i.e., `axis==0`; the rows) of `data` (and similarly for `label`) correspond 
    to the theta (polar) direction, and the second axis (i.e., `axis==1`; the columns)
    correspond to the azimuthal direction. The last column is simply the same as 
    the first (i.e., `data[0, :] == data[-1, :]`). The first and last rows must 
    represent the north and south poles such that each data[:, 0] and data[:, -1]
    slice contain the same values. Otherwise, unpredictable behavior may arise.
    """
    
    # ensuring the input is of the expected type for C interface
    data    = data.astype(dtype=_FLOAT, order='C', copy=False)
    labels  = labels.astype(dtype=_INT, order='C', copy=False)


    # total size of the problem
    N = grid.map.contents.num_nodes
    num_edges = grid.map.contents.num_edges

    # vectors to contain non-diagonal Laplacian entries
    rows    = np.zeros(num_edges*2, dtype=_INT,   order='C')
    columns = np.zeros(num_edges*2, dtype=_INT,   order='C')
    values  = np.zeros(num_edges*2, dtype=_FLOAT, order='C')

    log.debug("Calling get_ordered_Laplacian_from_map form C")
    num_labelled = platerecipy_clib_segmentation.get_ordered_Laplacian_from_map(
        grid.map,
        data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        ctypes.c_double(beta),
        rows.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        columns.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        values.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )
    
    L_ord = sp.sparse.csr_matrix((values, (rows, columns)), shape=(N, N), dtype=_FLOAT)
    L_ord.setdiag(-np.ravel(L_ord.sum(axis=0)))
    
    L_U = L_ord[num_labelled:, num_labelled:]
    B_T = L_ord[num_labelled:, :num_labelled]

    M = np.zeros((num_labelled, largest_label), dtype=_FLOAT, order='C')
    log.debug("Calling get_ordered_boundary_matrix_from_map form C")
    platerecipy_clib_segmentation.get_ordered_boundary_matrix_from_map(
        grid.map,
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        ctypes.c_int32(largest_label),
        ctypes.c_int32(num_labelled),
        M.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    
    A = L_U
    RHS = -B_T@M

    log.debug("The linear system is constructed.")
    return A, RHS, num_labelled

def solve_AX_B(
    A           : np.ndarray, 
    B           : np.ndarray,
    tol         = 1e-3,
    solver      = 'LU', 
    safe_mode   = True
) -> np.ndarray:
    """
    Solving for the sparse linear system arising from RW.

    Parameters
    ----------
    A : np.ndarray,
        The sparse left-hand-side matrix.

    B : np.ndarray,
        The right-hand-side matrix (different from B in RW) and must have at 
        least two columns.

    tol : float, default=1e-3,
        Tolerance values for `solver` options that require one (e.g., `solver='CG'`).
    
    solver : str, default='LU', 
        Choice of sparse solver. Possible options are: `'direct'`, `'LU'`, 
        `'CG'`, or `'FA'`.

    safe_mode : bool, default=True
        Whether to solver for the last column or indirectly solve for it noting 
        that the sum along each row must equate to one.
    
    
    Returns
    -------
    np.ndarray
    """
    
    n_RHS = B.shape[1]      # number of right-hand-sides
    X = np.zeros_like(B)    # same shape as B

    if not safe_mode:
        # not solving for the last label since the sum of probabilities should 
        # add up to zero
        log.warning("The safe mode is turned off. Now calculating the last columns (i.e., label) given the probabilities should add up to zero.")
        n_RHS -= 1


    if solver == 'direct':
        # solving Ax = b directly 
        for j in range(n_RHS): 
            X[:, j] = sp.sparse.linalg.spsolve(A, B[:, j])
    elif solver == 'CG':
        # solving it using conjugate gradient
        for j in range(n_RHS): 
            X[:, j], info = sp.sparse.linalg.cg(A, B[:, j], tol)
            if info > 0:
                log.warning(f'Failed to converge using {solver} solver with a tolerance of {tol} for column {j+1}.')
    elif solver == 'LU':
        # solving by LU decomposition
        lu = sp.sparse.linalg.splu(A)       # using the same LU for all for efficiency 
        for j in range(n_RHS): 
            X[:, j] = lu.solve(B[:, j])
    elif solver == 'FA':
        # solving using factorization
        FA_solve = sp.sparse.linalg.factorized(A)
        for j in range(n_RHS):
            X[:, j] = FA_solve(B[:, j])
    
    #elif solver == 'MG':
    #    ml = pyamg.smoothed_aggregation_solver(A)  # construct the multigrid hierarchy
    #    log.debug(str(ml))                         # print hierarchy information
    #    for j in range(n_RHS):
    #        X[:,j], info = ml.solve(B[:,j], tol=tol, return_info=True, cycle='V')                          # solve Ax=b to a tolerance of 1e-10
    #        if info > 0:
    #            log.warning(f'Failed to converge using {solver} solver with a tolerance of {tol} for column {j+1}.')
    
    else:
        raise ValueError('Unknown solver.')
    
    if safe_mode:
        X /= np.sum(X, axis=1)[:, None]
    else:
        X[:, -1] = 1. - np.sum(X[:, :-1], axis=1)

    return X


def random_walker(
    data                            : np.ndarray, 
    labels                          : np.ndarray,
    grid                            : Grid,
    beta                            = 100.,
    solver                          = 'direct',
    solver_tol                      = 1e-3,
    max_beta_reduction_attempts     = 2,
    beta_reduction_fraction         = 0.9,
    beta_reduction_residual_tol     = 1e-3
) -> tuple:
    """
    Perform Random Walker segmentation on a 2D input `data` field.

    Parameters
    ----------
    data : np.ndarray, 
        The input 2D gridded field.
    
    labels : np.ndarray,
        The labelled field. Must be the same shape as `data` and integer valued.
        The unmarked must be zeros and each marked patch must be identified with 
        a distinct positive integer. The marking integers must be sequential.

    beta : float, default=100.,
        The Gaussian scaling factor for the connection weights. This controls for 
        the sharpness of the identified segments, but with the trade off of 
        possible destabilization of the linear system.

    grid : Grid, optional
        An instance of the grid object to perform segmentation on.
     
    solver : str, default='direct',
        Choice of sparse linear solver. Possible options are:
        `'LU'`, `'FA'`, `'CG'`, and `'direct'`.

    solver_tol : float, default=1e-3,
        For the solvers that require a tolerance (e.g., if `solver='CG'`).

    max_beta_reduction_attempts : int, default=2,
        Maximum number of attempts to reduce `beta` in case the choice of `beta` 
        results in an ill-conditioned matrix or residual exceeds the value 
        specified by `beta_reduction_residual_tol`. To disable this mechanism, 
        this argument is to be set at 0.

    beta_reduction_fraction : float, default=0.9,
        The faction to reduce `beta` by (this value is multiplied to `beta` at 
        each attempt).

    beta_reduction_residual_tol : float, default=1e-3,
        The tolerance used to assessed whether the RW solution is converged.

        
    Returns
    -------
    tuple(`IDs`: np.ndarray, `probs`: np.ndarray)

    Warning
    -------
    If `is_spherical == True`, it is assumed that the first axis 
    (i.e., `axis==0`; the rows) of `data` (and similarly for `label`) correspond 
    to the theta (polar) direction, and the second axis (i.e., `axis==1`; the columns)
    correspond to the azimuthal direction. The last column is simply the same as 
    the first (i.e., `data[0, :] == data[-1, :]`). The first and last rows must 
    represent the north and south poles such that each data[:, 0] and data[:, -1]
    slice contain the same values. Otherwise, unpredictable behavior may arise.
    """
    # ensuring the input is of the expected type for C interface
    data    = data.astype(dtype=_FLOAT, order='C', copy=False)
    labels  = labels.astype(dtype=_INT, order='C', copy=False)

    unique_labels = np.unique(labels.ravel())
    if unique_labels.max()+1 == unique_labels.size:
        largest_label = unique_labels.max()
    else:
        log.error("Labels must be sequential positive integers starting from 1.")
        raise ValueError("Labels are not sequential.")

    if largest_label == 1:
        log.warning('Only 1 label was found.')
        return np.ones(data.shape, dtype=_INT, order='C'), \
            np.ones(data.shape, dtype=_FLOAT, order='C')
    
    # since large values of beta may destabilize the numerical system (such
    # that the solution gives rise to sum of probabilities greater than 1)
    # a built-in mechanism is implemented that systematically reduces beta
    # until convergence is reached.
    attempt = 0
    while attempt < max_beta_reduction_attempts + 1:
        log.debug(f"Constructing the linear system with beta={beta * (beta_reduction_fraction**attempt)}.")
        A, B, num_labelled = get_AB_from_map(
            data             = data,
            labels           = labels,
            beta             = beta * (beta_reduction_fraction**attempt),
            largest_label    = largest_label,
            grid             = grid
        )

        log.debug("Solving the linear system.")
        X = solve_AX_B(A, B, tol=solver_tol, solver=solver)

        residual = np.linalg.norm(B-A@X)
        log.info(f"Residual of the linear solution = {residual:4.2e}.")

        if residual < beta_reduction_residual_tol:
            break    

        attempt += 1

    log.debug("Generating IDs and probabilities from X.")
    
    X       = X.astype(dtype=_FLOAT, order='C', copy=False)
    IDs     = np.zeros(data.shape, dtype=_INT, order='C')
    probs   = np.zeros((*data.shape, X.shape[1]), _FLOAT, order='C')
    
    log.debug("Calling get_IDs_and_probs_from_X_and_map from C")
    platerecipy_clib_segmentation.get_IDs_and_probs_from_X_and_map(
        grid.map,
        X.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int32(X.shape[0]),
        ctypes.c_int32(X.shape[1]),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        ctypes.c_int32(num_labelled),
        IDs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        probs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ) 

    # constructing the composite probability field
    if len(probs.shape) == 2:
        # one-dimensional data
        i_max, _ = probs.shape
        ID_probs = np.zeros(i_max, dtype=_FLOAT, order='C')

        for i in range(i_max):
            prob_vec    = probs[i, :]
            ID_probs[i] = prob_vec[np.argmax(prob_vec)]
    
    elif len(probs.shape) == 3:
        # two-dimensional data
        i_max, j_max, _ = probs.shape
        ID_probs = np.zeros((i_max, j_max), dtype=_FLOAT, order='C')

        for i in range(i_max):
            for j in range(j_max):
                prob_vec    = probs[i, j, :]
                ID_probs[i, j] = prob_vec[np.argmax(prob_vec)]
    
    else:
        log.error("Unrecognized probability strucutre.")
        raise ValueError("Unrecognized probability strucutre.")

    if isinstance(grid, SphericalGrid):
        grid.enforce_data_consistency(IDs)
        grid.enforce_data_consistency(ID_probs)
    
    return IDs, ID_probs

