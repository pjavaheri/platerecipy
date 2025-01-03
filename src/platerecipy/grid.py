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

class Grid(object):
    """
    Generic class for a grid structure.
    """
    def __init__(self) -> Grid:
        '''Construct a generic grid structure.'''
        self._xs = None
        self._ys = None
        self._zs = None

    @property
    def xs(self):
        '''2D mesh grid of Cartesian x coordinates'''
        return self._xs

    #@xs.setter
    #def xs(self, value: np.ndarray):
    #    self._xs = value
    
    @property
    def ys(self):
        '''2D mesh grid of Cartesian y coordinates'''
        return self._ys
    
    @property
    def zs(self):
        '''2D mesh grid of Cartesian z coordinates'''
        return self._zs
    

class SphericalGrid(Grid):
    """
    Uniform spherical grid.
    """

    def __init__(
        self,
        original_xs  : np.ndarray,
        original_ys  : np.ndarray,
        original_zs  : np.ndarray,
        theta_res    = None,
        phi_res      = None,    
    ) -> SphericalGrid:
        """
        Constructs a uniform spherical grid structure. 
        
        If `theta_res` and `phi_res` are provided, the resolution will be 
        determined by those values. Alternatively, the resolution close enough 
        to what provided by the Cartesian arguments.

        Warning
        -------
        The Cartesian arguments will be the basis upon which the class method
        `interpolate_field()` operates. In consistent values will impact the 
        interpolation process.

        Parameters
        ----------
        original_xs : np.ndarray,
            array of original input Cartesian x coordinates.

        original_ys : np.ndarray,
            array of original input Cartesian y coordinates.

        original_zs : np.ndarray,
            array of original input Cartesian z coordinates.
        
        theta_res : int, optional
            polar (colatitude) resolution. 

        phi_res : int, optional
            azimuthal (longitude) resolution.
        """
        Grid.__init__(self)
        self._r = None
        ## the order is colatitude (theta) and azimuth (phi)

        self.original_xs = original_xs
        self.original_ys = original_ys
        self.original_zs = original_zs

        if phi_res is None:
            # number of grid points along the azimuth
            phi_res = int(np.sqrt(2 * self.original_xs.size))
        
        if theta_res is None:
            # number of grid points along the colatitude
            theta_res = int(np.sqrt(self.original_xs.size / 2))
        
        self.theta_res = theta_res
        self.phi_res   = phi_res

        self.original_rs     = np.sqrt(
            self.original_xs**2 + self.original_ys**2 + self.original_zs**2
        )
        self.original_thetas = np.arccos(self.original_zs/self.original_rs)
        self.original_phis   = np.arctan2(self.original_ys, self.original_xs)
        
        self.original_points_in_theta_phi = np.vstack(
            [
                self.original_thetas.ravel(), self.original_phis.ravel() # flatten instead?
            ]
        ).T

        self._mercator_thetas = np.linspace(0, np.pi, self.theta_res)
        self._mercator_phis = np.linspace(-np.pi, np.pi, self.phi_res)
        self._thetas, self._phis = np.meshgrid(
            self._mercator_thetas, 
            self._mercator_phis, # there is an issue here. incorrect phi
            indexing='ij'
        )

        self._r = self.original_rs.mean()

        self._xs = self._r*np.sin(self._thetas)*np.cos(self._phis)
        self._ys = self._r*np.sin(self._thetas)*np.sin(self._phis)
        self._zs = self._r*np.cos(self._thetas)

    @property
    def r(self):
        '''Radius at the surface'''
        return self._r
    
    @property
    def thetas(self):
        '''2D grid of spherical colatitude coordinates [0, pi]'''
        return self._thetas
    
    @property
    def phis(self):
        '''2D grid of spherical azimuthal coordinates [-pi, pi]'''
        return self._phis

    def __str__(self):
        return f'A Grid object with\ninterp_theta_res={self.interp_theta_res}\ninterp_phi_res={self.interp_phi_res}\n'


    def interpolate_field(
        self,
        field       : np.ndarray,
        take_log    = False
    ) -> np.ndarray:
        """
        Interpolates the input field linearly for convex hull interior points and
        closest points to the exterior.

        Parameters
        ----------
        field : np.ndarray,
            array corresponding to the original Cartesian coordinate inputs.
        
        take_log : bool, False
            whether to perform the interpolation on the log space (suitable for
            fields that vary by orders of magnitude), default = False.
        
        Returns
        -------
        np.ndarray
        """
        field = field.ravel()

        # checking if interpolation needs to occur on the log space
        if take_log:
            field = np.log(field)

        # this will leave the boundary (out of the convex hull) as NaNs
        grid_linear = griddata(
            self.original_points_in_theta_phi, 
            field, 
            (self.thetas, self.phis), 
            method='linear'
        )

        # replacing NaNs with closest values
        grid_nearest = griddata(
            self.original_points_in_theta_phi, 
            field, 
            (self.thetas, self.phis), 
            method='nearest'
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
        
        return gridded_field
