"""
@file output.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for grid input, output and visualizations.
"""
import logging
log = logging.getLogger(__name__)

import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import pyvista as pv
from pandas import DataFrame

from . import _FLOAT

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
    log.debug("Converting the grid to a pyvsita mesh ...")
    # setting up the mesh according to the first field
    gridded_field, field_name = gridded_fields[0], field_names[0]
    # since pyvista wants [0,360] not [-90,90]
    gridded_field = np.roll(gridded_field, gridded_field.shape[1]//2, axis=1)

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
            # since pyvista wants [0,360] not [-90,90]
            gridded_field = np.roll(gridded_field, gridded_field.shape[1]//2, axis=1)
 
            north_pole = gridded_field[ 0, :].mean(dtype=gridded_field.dtype)
            south_pole = gridded_field[-1, :].mean(dtype=gridded_field.dtype)
            gridded_field = gridded_field[1:-1, :].flatten('F')
            
            gridded_field = np.hstack([north_pole, south_pole, gridded_field])
            
            # adding the field
            sph_mesh[field_name] = gridded_field
    
    log.debug('... mesh created.')
    
    return sph_mesh

class Grid(object):
    """
    Generic class for a grid structure.
    """
    def __init__(self):
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
    

class PartialSphericalGrid(Grid):
    """
    Partial uniform spherical grid.
    """

    def __init__(
        self,
        original_xs  : np.ndarray,
        original_ys  : np.ndarray,
        original_zs  : np.ndarray,
        theta_range  : tuple,
        phi_range    : tuple,
        theta_res    = None,
        phi_res      = None,    
    ):
        """
        Constructs a partial uniform spherical grid structure, spanning  
        [`theta_range[0]`,`theta_range[1]`] x [`phi_range[0]`, `phi_range[1]`].
        
        If `theta_res` and `phi_res` are provided, the resolution will be 
        determined by those values. Alternatively, the resolution close enough 
        to what provided by the Cartesian arguments.

        Warning
        -------
        The Cartesian arguments will be the basis upon which the class method
        `interpolate_field()` operates. In consistent values will impact the 
        interpolation process.

        Warning
        -------
        The input requires a range that is increasing in phi. That is, ranges 
        such that cross the -pi/+pi boundary are not acceptable unless the full
        sphere is covered. An examples of unacceptable range is [pi/2, -pi/3].
        Max range for theta: [0, pi]
        Max range for phi: [-pi, pi]

        Parameters
        ----------
        original_xs : np.ndarray,
            array of original input Cartesian x coordinates.

        original_ys : np.ndarray,
            array of original input Cartesian y coordinates.

        original_zs : np.ndarray,
            array of original input Cartesian z coordinates.
        
        theta_range : tuple,
            the polar (colatitude) range.
        
        phi_range : tuple,
            the azimuthal (longitude) range.
        
        theta_res : int, optional
            polar (colatitude) resolution. 

        phi_res : int, optional
            azimuthal (longitude) resolution.
        """
        log.debug('Generating a PartialSphericalGrid ...')

        Grid.__init__(self)
        self._r = None
        ## the order is colatitude (theta) and azimuth (phi)

        self.original_xs = original_xs.astype(dtype=_FLOAT, order='C', copy=False)
        self.original_ys = original_ys.astype(dtype=_FLOAT, order='C', copy=False)
        self.original_zs = original_zs.astype(dtype=_FLOAT, order='C', copy=False)

        self._theta_range = theta_range
        self._phi_range   = phi_range

        range_ratio = (phi_range[1]-phi_range[0])/(theta_range[1]-theta_range[0])

        if phi_res is None:
            # number of grid points along the azimuth
            phi_res = int(np.sqrt(range_ratio * self.original_xs.size))
        
        if theta_res is None:
            # number of grid points along the colatitude
            theta_res = int(np.sqrt(self.original_xs.size / range_ratio))
        
        self._theta_res = theta_res
        self._phi_res   = phi_res

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

        self._mercator_thetas = np.linspace(
            self._theta_range[0], 
            self._theta_range[1], 
            self._theta_res, dtype=_FLOAT
        )
        self._mercator_phis = np.linspace(
            self._phi_range[0], 
            self._phi_range[1], 
            self._phi_res, dtype=_FLOAT
        )
        self._thetas, self._phis = np.meshgrid(
            self._mercator_thetas, 
            self._mercator_phis, # there is an issue here. incorrect phi
            indexing='ij'
        )

        self._r = self.original_rs.mean()

        self._xs = self._r*np.sin(self._thetas)*np.cos(self._phis)
        self._ys = self._r*np.sin(self._thetas)*np.sin(self._phis)
        self._zs = self._r*np.cos(self._thetas)

        log.debug('... grid created.')

    @property
    def r(self):
        '''Radius at the surface'''
        return self._r

    @property
    def theta_range(self):
        '''Range of theta (polar angle)'''
        return self._theta_range

    @property
    def phi_range(self):
        '''Range of phi (polar angle)'''
        return self._phi_range
    
    @property
    def thetas(self):
        '''2D grid of spherical colatitude coordinates [0, pi]'''
        return self._thetas
    
    @property
    def phis(self):
        '''2D grid of spherical azimuthal coordinates [-pi, pi]'''
        return self._phis

    def __str__(self):
        return f'A PartialSphericalGrid object with\ntheta_range={self.theta_range}\nphi_range={self.phi_range}\ninterp_theta_res={self.interp_theta_res}\ninterp_phi_res={self.interp_phi_res}\n'


    def interpolate_field(
        self,
        field           : np.ndarray,
        take_log        = False,
        method          = "lat-lon", #"lat-lon" or "tangent-plane",
        **method_kwargs
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
        log.debug('Interpolating the data to a uniform spherical grid ...')


        field = field.ravel().astype(dtype=_FLOAT, order='C', copy=False)

        # checking if interpolation needs to occur on the log space
        if take_log:
            field = np.log(field)

        if method == "lat-lon":
            log.debug('... lat-lon linear interpolation')

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

        elif method == "tangent-plane":
            log.debug('... tangent linear interpolation using an LSQ fit to the closest neighbors')

            # loading defaults
            if not 'k' in method_kwargs.keys():
                method_kwargs['k'] = 4
            
            if not 'eps' in method_kwargs.keys():
                method_kwargs['eps'] = 1e-6
            
            if not 'near_thresh' in method_kwargs.keys():
                method_kwargs['near_thresh'] = 1e-6
            
            if not hasattr(self, '_original_cart_mat'):
                self._original_cart_mat = np.vstack([
                    self.original_xs.ravel(), 
                    self.original_ys.ravel(), 
                    self.original_zs.ravel(),
                ]).T
            
            if not hasattr(self, '_cart_mat'):
                self._cart_mat = np.vstack([
                    self.xs.ravel(), 
                    self.ys.ravel(), 
                    self.zs.ravel(),
                ]).T
            
            # constructing a similarity tree
            if not hasattr(self, '_original_kdtree'):
                log.debug('... constructing a similarity tree')
                self._original_kdtree = cKDTree(self._original_cart_mat)
            
            if (not hasattr(self, '_original_neigh_dists')) or (not hasattr(self, '_original_neighs')):
                log.debug('... querying the tree for near neighbors')
                self._original_neigh_dists, self._original_neighs = \
                    self._original_kdtree.query(self._cart_mat, k=method_kwargs['k'])

            neighs     = self._original_cart_mat[self._original_neighs]
            neigh_vals = field[self._original_neighs]            

            # unit vectors on the surface
            #r_hats = self._original_cart_mat / self.r
            theta_hats = np.vstack([
                np.cos(self.thetas.ravel())*np.cos(self.phis.ravel()),
                np.cos(self.thetas.ravel())*np.sin(self.phis.ravel()),
                -np.sin(self.thetas.ravel())
            ]).T
            phi_hats = np.vstack([
                -np.sin(self.phis.ravel()),
                np.cos(self.phis.ravel()),
                np.zeros_like(self.thetas.ravel())
            ]).T

            # constructing a linear model on the tangent space
            log.debug('... constructing a linear LSQ model on the tangent space')
            ps = np.array([
                neighs[:, i, :] - self._cart_mat[:, :] \
                    for i in range(method_kwargs['k'])
            ])
            fs = [neigh_vals[:, i] for i in range(method_kwargs['k'])]
            xs = np.array([np.vecdot(ps[i], theta_hats, axis=1) for i in range(method_kwargs['k'])])
            ys = np.array([np.vecdot(ps[i], phi_hats, axis=1) for i in range(method_kwargs['k'])])
            G = np.array([
                [np.ones_like(xs[i]), xs[i], ys[i]] for i in range(method_kwargs['k'])
            ])
            G = np.transpose(G, [2, 0, 1])
            GT = np.transpose(G, [0, 2, 1])
            GTG = GT@G
            trace_GTG = np.trace(GTG, axis1=1, axis2=2)
            GTG[:, 0, 0] += method_kwargs['eps']*trace_GTG
            GTG[:, 1, 1] += method_kwargs['eps']*trace_GTG
            GTG[:, 2, 2] += method_kwargs['eps']*trace_GTG
            B = np.vstack([fs[i] for i in range(method_kwargs['k'])])
            B = np.transpose(B, [1, 0])
            B = B.reshape((B.shape[0], B.shape[1], 1))

            log.debug('... solving all linear LSQ problems together (stacked)')
            grid_linear = (np.linalg.solve(GTG, GT@B)[:, 0]).ravel() # this may fail

            # if a point is too close to a gride point, simply use its value
            grid_linear[
                self._original_neigh_dists[:, 0]/self.r < method_kwargs['near_thresh']
            ] = fs[0][
                self._original_neigh_dists[:, 0]/self.r < method_kwargs['near_thresh']
            ]

            grid_linear = grid_linear.reshape(self.xs.shape)
        
        # forcing boundary conditions by averaging
        if np.abs(self.theta_range[0]) < 1e-8:
            grid_linear[ 0,  :] = np.mean(grid_linear[:,  0])
        if np.abs(self.theta_range[1] - np.pi) < 1e-8:
            grid_linear[-1,  :] = np.mean(grid_linear[:, -1])
        if (np.abs(-np.pi - self.phi_range[0]) < 1e-8) \
            and (np.abs(np.pi - self.phi_range[1]) < 1e-8):
            grid_linear[ :,  0] = 0.5 * (grid_linear[:, 0] + grid_linear[:, -1])
            grid_linear[ :, -1] = grid_linear[:, 0]
        
        # returning to original space
        if take_log:
            gridded_field = np.exp(grid_linear).astype(dtype=_FLOAT, order='C', copy=False)
        else:
            gridded_field = grid_linear.astype(dtype=_FLOAT, order='C', copy=False)
        
        log.debug('... interpolation done.')
        return gridded_field

    def map_to_original_input(
        self,
        field       : np.ndarray,
        method      = 'nearest',
        take_log    = False,
        csv_output  = None
    ):
        """
        Maps the spherical field back to the original input format by interpolating 
        the nodes using the nearest method.

        Parameters
        ----------
        field : np.ndarray,
            array corresponding to the original Cartesian coordinate inputs.
        
        method : str, optional,
            'nearest' or 'linear', default = 'nearest'

        take_log : bool, False,
            whether to perform the interpolation on the log space (suitable for
            fields that vary by orders of magnitude), and makes a substantive
            difference if `method == 'linear'`, default = False.
        
        csv_output : str, optional,
            whether to save the data in a CSV file with name csv_output
        
        Returns
        -------
        np.ndarray
        """
        # checking if interpolation needs to occur on the log space
        if take_log:
            field = np.log(field)

        org_nearest = griddata(
            (self.thetas.ravel(), self.phis.ravel()),
            field.ravel(),  
            self.original_points_in_theta_phi, 
            method='nearest'
        )

        if method == 'linear':
            org_linear = griddata(
                (self.thetas.ravel(), self.phis.ravel()),
                field.ravel(),  
                self.original_points_in_theta_phi, 
                method='linear'
            )
            org_linear[np.isnan(org_linear)] = org_nearest[np.isnan(org_linear)]
            org_order = org_linear
        else:
            org_order = org_nearest

        if take_log:
            org_order = np.exp(org_order)
        
        if csv_output is not None:
            data = {
                'theta' : self.original_points_in_theta_phi[:, 0],
                'phi'   : self.original_points_in_theta_phi[:, 1],
                'value' : org_order 
            }
            DataFrame(data).to_csv(csv_output)

        return org_order



class SphericalGrid(PartialSphericalGrid):
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
    ):
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
        log.debug('Generating a SphericalGrid ...')

        PartialSphericalGrid.__init__(
            self, 
            original_xs = original_xs,
            original_ys = original_ys,
            original_zs = original_zs,
            theta_range = [0., np.pi],
            phi_range   = [-np.pi, np.pi],
            theta_res   = theta_res,
            phi_res     = phi_res
        )

    def __str__(self):
        return f'A SphericalGrid object with\ninterp_theta_res={self.interp_theta_res}\ninterp_phi_res={self.interp_phi_res}\n'
