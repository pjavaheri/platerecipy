"""
@file output.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for grid input, output and visualizations.
"""

import numpy as np
from scipy.interpolate import griddata
import pyvista as pv

def interpolate_to_spherical(
    xs          : np.ndarray,
    ys          : np.ndarray,
    zs          : np.ndarray,
    values      : np.ndarray,
    phi_res     = None,
    theta_res   = None,
    take_log    = False,
    is_scalar   = True,
) -> tuple:
    ## the order is colatitude (theta) and azimuth (phi)

    if phi_res is None:
        # number of grid points along the azimuth
        phi_res = int(np.sqrt(2 * xs.size))
    
    if theta_res is None:
        # number of grid points along the colatitude
        theta_res = int(np.sqrt(xs.size / 2))
    
    
    xs = xs.ravel()
    ys = ys.ravel()
    zs = zs.ravel()

    if is_scalar:
        values = values.ravel()
    else:
        values = values.reshape(xs.size, 3)

    rs      = np.sqrt(xs**2 + ys**2 + zs**2)
    thetas  = np.arccos(zs/rs)
    phis    = np.arctan2(ys, xs)
    
    points  = np.vstack([phis, thetas]).T
    
    grid_theta, grid_phi = np.meshgrid(
        np.linspace(0, np.pi, theta_res), 
        np.linspace(-np.pi, np.pi, phi_res), # there is an issue here. incorrect phi
        indexing='ij'
    )

    # performing linear interpolation

    # checking if interpolation needs to occur on the log space
    if take_log:
        values = np.log(values)

    # this will leave the boundary (out of the convex hull) as NaNs
    grid_linear = griddata(
        points, values, (grid_phi, grid_theta), method='linear'
    )

    # replacing NaNs with closest values
    grid_nearest = griddata(
        points, values, (grid_phi, grid_theta), method='nearest'
    )
    grid_linear[np.isnan(grid_linear)] = grid_nearest[np.isnan(grid_linear)]
    
    # forcing boundary conditions by averaging
    grid_linear[ 0,  :] = np.mean(grid_linear[:,  0])
    grid_linear[-1,  :] = np.mean(grid_linear[:, -1])
    grid_linear[ :,  0] = 0.5 * (grid_linear[:, 0] + grid_linear[:, -1])
    grid_linear[ :, -1] = grid_linear[:, 0]
    
    # returning to original space
    if take_log:
        gridded_field = np.exp(grid_linear)
    else:
        gridded_field = grid_linear
    
    return gridded_field, (rs.mean(), grid_theta, grid_phi)

def convert_grid_to_mesh(
    gridded_fields  : list,
    field_names     : list,
    radius          = 1.
):
    """
    Stores all gridded fields onto a spherical mesh (with nodal connectivity).
    
    Parameters
    ----------
    gridded_fields  : list,
    field_names     : list,
    theta_phi_axes  = None,
    radius          = 1.

    Returns
    -------
    pyvista mesh

    Warning
    -------
    Gridded fields are assumed to be C-like (row-major) arrays with the first 
    dimension corresponding to theta (the polar angle) and the second to phi 
    (the azimuthal angle).

    """
    # setting up the mesh according to the first field
    gridded_field, field_name = gridded_fields[0], field_names[0]

    theta_resolution  = gridded_field.shape[0]
    phi_resolution    = gridded_field.shape[1]

    # unravelling the 2D field into 1D and collapsing the polar singularities
    north_pole = gridded_field[ 0, :].mean(dtype=gridded_field.dtype)
    south_pole = gridded_field[-1, :].mean(dtype=gridded_field.dtype)
    gridded_field = gridded_field[1:-1, :].flatten('F')

    # adding singular poles to the gridded field
    gridded_field = np.hstack([north_pole, south_pole, gridded_field])
        
    # pyvista's definition for theta and phi is opposite
    sph_mesh = pv.Sphere(
        radius              = radius, 
        theta_resolution    = phi_resolution, 
        phi_resolution      = theta_resolution
    )

    # adding the field
    sph_mesh[field_name] = gridded_field

    # if more than one field
    if len(gridded_fields) > 1:
        for i in range(1, len(gridded_fields)):
            gridded_field, field_name = gridded_fields[i], field_names[i]
 
            north_pole = gridded_field[ 0, :].mean(dtype=gridded_field.dtype)
            south_pole = gridded_field[-1, :].mean(dtype=gridded_field.dtype)
            gridded_field = gridded_field[1:-1, :].flatten('F')
            
            gridded_field = np.hstack([north_pole, south_pole, gridded_field])
            
            # adding the field
            sph_mesh[field_name] = gridded_field
    
    return sph_mesh


