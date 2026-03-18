"""
@file model.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for model post-processing and analysis.
"""

import numpy as np
from .grid import SphericalGrid
from .model import PlateModel
from .transform import single_plate_interior_distance_transform, \
    full_plate_interior_distance_transform
from . import _INT, _FLOAT

# only spherical distances for plateness
# planar needs to be defined

def get_plateness(
    edots       : np.ndarray,
    model       : PlateModel,
    f           = 0.8,
    plate_ID    = None,
    num_threads = 1  # or 'auto'
):
    edots = edots.astype(dtype=_FLOAT, order='C', copy=False)
    
    plateness = {}
    dA = (model.grid.r**2) * np.sin(model.grid.thetas)
    log_edots = np.log(edots)  

    if model.plate_IDs is None:
        raise ValueError("PlateModel object must be first called to find its plates before passing for analysis.")
    
    if plate_ID is not None:
        # calculating plateness for a single plate ID
        plate_mask = model.plate_IDs == plate_ID

        if isinstance(model.grid, SphericalGrid):
            # a spherical distance transform
            D = single_plate_interior_distance_transform(
                xs=model.grid.xs,
                ys=model.grid.ys,
                zs=model.grid.zs,
                plate_indicators=model.plate_IDs == plate_ID,
                R=model.grid.r,
                num_threads=num_threads
            )          
            
            edots_thresh_local = \
                np.min(edots[plate_mask]) \
                    + f*(np.max(edots[plate_mask]) - np.min(edots[plate_mask]))
            
            log_edots_thresh_global = \
                log_edots.min() + f*(log_edots.max() - log_edots.min())
            log_edots_thresh_local = \
                np.min(log_edots[plate_mask]) \
                    + f*(np.max(log_edots[plate_mask]) - np.min(log_edots[plate_mask]))
 
            mask = (plate_mask) & (edots > edots_thresh_local)
            plateness['legacy_local'] = \
                1. - np.sum(dA[mask]) / np.sum(dA[plate_mask]) / 0.6
            
            mask = (plate_mask) & (log_edots > log_edots_thresh_local)
            plateness['modified_local'] = \
                1. - np.sum(D[mask]*dA[mask]) / np.sum(D[plate_mask]*dA[plate_mask])

            mask = (plate_mask) & (log_edots > log_edots_thresh_global)
            plateness['modified_global'] = \
                1. - np.sum(D[mask]*dA[mask]) / np.sum(D[plate_mask]*dA[plate_mask])

        else:
            # planar
            pass

    else:
        # calculating plateness for all identified plates
        unique_IDs = np.unique(model.plate_IDs)
        num_plates = len(unique_IDs)

        plateness['legacy_local']      = np.zeros(num_plates, _FLOAT, order='C')
        plateness['modified_local']    = np.zeros(num_plates, _FLOAT, order='C')
        plateness['modified_global']   = np.zeros(num_plates, _FLOAT, order='C')

        if isinstance(model.grid, SphericalGrid):
            # a spherical distance transform
            D = full_plate_interior_distance_transform(
                xs=model.grid.xs,
                ys=model.grid.ys,
                zs=model.grid.zs,
                plate_IDs=model.plate_IDs,
                R=model.grid.r,
                num_threads=num_threads
            )

            log_edots_thresh_global = \
                log_edots.min() + f*(log_edots.max() - log_edots.min())

            for i, plate_ID in enumerate(unique_IDs):
                plate_mask = model.plate_IDs == plate_ID

                edots_thresh_local = \
                    np.min(edots[plate_mask]) \
                        + f*(np.max(edots[plate_mask]) - np.min(edots[plate_mask]))
                
                log_edots_thresh_local = \
                    np.min(log_edots[plate_mask]) \
                        + f*(np.max(log_edots[plate_mask]) - np.min(log_edots[plate_mask]))
    
                mask = (plate_mask) & (edots > edots_thresh_local)
                plateness['legacy_local'][i] = \
                    1. - np.sum(dA[mask]) / np.sum(dA[plate_mask]) / 0.6
                
                mask = (plate_mask) & (log_edots > log_edots_thresh_local)
                plateness['modified_local'][i] = \
                    1. - np.sum(D[mask]*dA[mask]) / np.sum(D[plate_mask]*dA[plate_mask])

                mask = (plate_mask) & (log_edots > log_edots_thresh_global)
                plateness['modified_global'][i] = \
                    1. - np.sum(D[mask]*dA[mask]) / np.sum(D[plate_mask]*dA[plate_mask])


        else:
            # planar
            pass
    
    return plateness

def get_euler_pole(
    vxs         : np.ndarray,
    vys         : np.ndarray,
    vzs         : np.ndarray,
    model       : PlateModel,
    plate_ID    = None,
) -> np.ndarray:
    if model.plate_IDs is None:
        raise ValueError("PlateModel object must be first called to find its plates before passing for analysis.")
    
    if plate_ID is not None:
        # calculating Euler poles for a single plate ID
        Px = model.grid.xs[model.plate_IDs == plate_ID].ravel()
        Py = model.grid.ys[model.plate_IDs == plate_ID].ravel()
        Pz = model.grid.zs[model.plate_IDs == plate_ID].ravel()
        PP = np.vstack([Px, Py, Pz]).T
        P = np.zeros((3*PP.shape[0], 3), dtype=_FLOAT)
        for i in range(PP.shape[0]):
            x, y, z = PP[i,:]
            P[3*i:3*(i+1), :] = np.array([
                [0, z, -y],
                [-z, 0, x],
                [y, -x, 0]
            ])[:,:]
        Vx = vxs[model.plate_IDs == plate_ID].ravel()
        Vy = vys[model.plate_IDs == plate_ID].ravel()
        Vz = vzs[model.plate_IDs == plate_ID].ravel()
        V = np.vstack([Vx, Vy, Vz]).T.ravel()
        R = np.linalg.solve(P.T@P, P.T@V).astype(dtype=_FLOAT, order='C', copy=False)
    else:
        # finding them for all plate IDs
        IDs = np.unique(model.plate_IDs)
        R = np.zeros((IDs.size, 3), dtype=_FLOAT, order='C')

        for i, plate_ID in enumerate(IDs):
            # calculating Euler poles for a single plate ID
            Px = model.grid.xs[model.plate_IDs == plate_ID].ravel()
            Py = model.grid.ys[model.plate_IDs == plate_ID].ravel()
            Pz = model.grid.zs[model.plate_IDs == plate_ID].ravel()
            PP = np.vstack([Px, Py, Pz]).T
            P = np.zeros((3*PP.shape[0], 3), dtype=_FLOAT)
            for i in range(PP.shape[0]):
                x, y, z = PP[i,:]
                P[3*i:3*(i+1), :] = np.array([
                    [0, z, -y],
                    [-z, 0, x],
                    [y, -x, 0]
                ])[:,:]
            Vx = vxs[model.plate_IDs == plate_ID].ravel()
            Vy = vys[model.plate_IDs == plate_ID].ravel()
            Vz = vzs[model.plate_IDs == plate_ID].ravel()
            V = np.vstack([Vx, Vy, Vz]).T.ravel()
            R[i, :] = np.linalg.solve(P.T@P, P.T@V).astype(dtype=_FLOAT, order='C', copy=False)
    
    return R

'''
V_rigid = np.zeros(PP.shape[0])
d = np.zeros(PP.shape[0])
for i in range(PP.shape[0]):
    V_rigid[i] = np.linalg.norm(np.cross(PP[i, :], R)[:])
    d[i] = np.linalg.norm(np.cross(PP[i, :], R)[:])/np.linalg.norm(R)
'''