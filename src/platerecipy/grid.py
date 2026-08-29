"""
File brief
----------
`grid.py`

Module for grid input, output and visualizations.

This is a part of `platerecipy` package. For license and citation, please
refer to the main repository:
[github.com/pjavaheri/platerecipy](github.com/pjavaheri/platerecipy)

Author(s): 
Pejvak Javaheri; [pejvak.javaheri@mail.utoronto.ca](mailto:pejvak.javaheri@mail.utoronto.ca)
"""

import logging
log = logging.getLogger(__name__)

import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from pathlib import Path

import pyvista as pv

from . import _FLOAT, _INT

# ~~~~~~~~~~~~~~ importing shared libraries ~~~~~~~~~~~~~
import weakref
import os
import sysconfig
import ctypes

module_path        : str  = os.path.dirname(os.path.abspath(__file__))
module_path        : str  = os.path.abspath(os.path.join(module_path, os.pardir))
shared_object_path : str  = os.path.join(
    module_path, 
    "libplaterecipy_grid" + sysconfig.get_config_var('EXT_SUFFIX')
)

log.debug("Loading platerecipy_clib_grid shared library")
"""
A Python access to `clib/grid.h` module.
"""
platerecipy_clib_grid = ctypes.CDLL(shared_object_path)

# ~~~~~~~~~~~~~~~ Grid type identifiers ~~~~~~~~~~~~~~~

# Custom/unknown grid
GTYPE_CUSTOM               = -1

# Flat rectilinear grid
GTYPE_PLANAR               = 0

# Partial rectilinear spherical grid
GTYPE_PARTIAL_SPHERICAL    = 1

# Spherical grid (with wraparound and polar nodes)
GTYPE_SPHERICAL            = 2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class Node(ctypes.Structure):
    _fields_ = [
        ("idx", ctypes.c_int32),
        ("npy_idx", ctypes.c_int32),
        ("ord_idx", ctypes.c_int32),
        ("neighs", ctypes.POINTER(ctypes.c_int32)),
        ("edge_lengths", ctypes.POINTER(ctypes.c_double)),
        ("num_neighs", ctypes.c_int32),
    ]


class Map(ctypes.Structure):
    _fields_ = [
        ("nodes", ctypes.POINTER(Node)),
        ("num_nodes", ctypes.c_int32),
        ("num_edges", ctypes.c_int32),
        ("i_max", ctypes.c_int32),
        ("j_max", ctypes.c_int32),
    ]


# platerecipy_clib_grid.get_node_at() function signature
platerecipy_clib_grid.get_node_at.argtypes = [
    ctypes.POINTER(Map),
    ctypes.c_int32
]
platerecipy_clib_grid.get_node_at.restype = ctypes.POINTER(Node)


# platerecipy_clib_grid.get_neigh_at() function signature
platerecipy_clib_grid.get_neigh_at.argtypes = [
    ctypes.POINTER(Node),
    ctypes.c_int32
]
platerecipy_clib_grid.get_neigh_at.restype = ctypes.c_int32


# platerecipy_clib_grid.set_npy_idx() function signature
platerecipy_clib_grid.set_npy_idx.argtypes = [
    ctypes.POINTER(Map),
    ctypes.POINTER(ctypes.c_int32)
]
platerecipy_clib_grid.set_npy_idx.restype = None


# platerecipy_clib_grid.alloc_map_from_cells() function signature
platerecipy_clib_grid.alloc_map_from_cells.argtypes = [
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_int32
]
platerecipy_clib_grid.alloc_map_from_cells.restype = ctypes.POINTER(Map)


# platerecipy_clib_grid.alloc_rect_map() function signature
platerecipy_clib_grid.alloc_rect_map.argtypes = [
    ctypes.c_int32, 
    ctypes.c_int32,
    ctypes.c_int32
]
platerecipy_clib_grid.alloc_rect_map.restype = ctypes.POINTER(Map)


# platerecipy_clib_grid.alloc_sph_map() function signature
platerecipy_clib_grid.alloc_sph_map.argtypes = [
    ctypes.c_int32, 
    ctypes.c_int32
]
platerecipy_clib_grid.alloc_sph_map.restype = ctypes.POINTER(Map)


# platerecipy_clib_grid.free_map() function signature
platerecipy_clib_grid.free_map.argtypes = [
    ctypes.POINTER(Map)
]
platerecipy_clib_grid.free_map.restype = None


# platerecipy_clib_grid.gen_cells_from_map() function signature
platerecipy_clib_grid.gen_cells_from_map.argtypes = [
    ctypes.POINTER(Map),
    ctypes.POINTER(ctypes.c_int32)
]
platerecipy_clib_grid.gen_cells_from_map.restype = None



# establishing a callback mechanism for enabling logging from c
callback_func_type = ctypes.CFUNCTYPE(None, ctypes.c_char_p)

def _grid_h_c_log(msg):
    log.debug(f"[platerecipy_clib_grid callback]: {msg.decode()}")

_grid_h_callback = callback_func_type(_grid_h_c_log) # not to be removed
platerecipy_clib_grid.set_grid_h_logger(_grid_h_callback)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



class Grid(object):
    """
    Generic class for a grid structure.
    """
    def __init__(self) -> None:
        '''Construct a generic grid structure.'''
        self._map   = None
        self._mesh  = None
        self._r     = None
        self._xs    = None
        self._ys    = None
        self._zs    = None
        
        self._finalizer = weakref.finalize(
            self,
            Grid._free_map,
            self._map
        )
    
    def enforce_data_consistency(self, array: np.ndarray) -> None:
        pass
        

    @property
    def map(self):
        '''Connection map of the grid'''
        return self._map

    @property
    def mesh(self):
        '''Pyvista mesh object for cross-platform io compatibility'''
        return self._mesh

    @property
    def r(self):
        '''Radius (or non-dimensional scaling factor)'''
        return self._r

    @r.setter
    def r(self, value: float):
        self._r = value

    @property
    def xs(self):
        '''2D mesh grid of Cartesian x coordinates'''
        return self._xs
    
    @property
    def ys(self):
        '''2D mesh grid of Cartesian y coordinates'''
        return self._ys
    
    @property
    def zs(self):
        '''2D mesh grid of Cartesian z coordinates'''
        return self._zs

    # deconstructor functions:
    @staticmethod
    def _free_map(map_ptr):
        if map_ptr is not None:
            log.debug("Deallocating Map struct")
            platerecipy_clib_grid.free_map(map_ptr)
    
    def __exit__(self, exc_type, exc, tb) -> None:
        self._finalizer()

    def __del__(self) -> None:
        self._finalizer()    


class CustomGrid(Grid):
    """
    Custom grid for the full sphere, defined by an input vtu file.
    """
    def __init__(
        self,
        surface_vtu_adr: str | Path,
    ):
        log.debug('Generating a CustomGrid from input ...')
        Grid.__init__(self)
        vtu_path = Path(surface_vtu_adr)

        if not vtu_path.is_file():
            raise FileNotFoundError(f'{vtu_path} not found.')
        
        mesh: pv.UnstructuredGrid = pv.read(vtu_path)

        self._mesh = mesh
        all_cells = mesh.cells.astype(_INT)
        
        self._map = platerecipy_clib_grid.alloc_map_from_cells(
            ctypes.c_int32(mesh.number_of_points),
            all_cells.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            ctypes.c_int32(mesh.cells.size)
        )

        self._xs = mesh.points[:, 0].astype(_FLOAT)
        self._ys = mesh.points[:, 1].astype(_FLOAT)
        self._zs = mesh.points[:, 2].astype(_FLOAT)

        # assuming that it's spherical; otherwise it would be used for normalizing distances
        self._r = np.sqrt(self._xs**2 + self._ys**2 + self._zs**2).mean()

        # assigning indices 
        npy_idxs = np.linspace(0, mesh.number_of_points-1, mesh.number_of_points, dtype=_INT)
        log.debug("Calling set_npy_idx from C")
        platerecipy_clib_grid.set_npy_idx(
            self._map,
            npy_idxs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        )

        log.debug("Calling set_cartesian_edge_length from C")
        platerecipy_clib_grid.set_cartesian_edge_length(
            self._map,
            self._xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_double(self._r)
        )

        log.debug('... grid created.')




class PartialSphericalGrid(Grid):
    """
    Partial uniform spherical grid.
    """

    def __init__(
        self,
        original_xs  : np.ndarray,
        original_ys  : np.ndarray,
        original_zs  : np.ndarray,
        theta_range  : tuple[float, float],
        phi_range    : tuple[float, float],
        theta_res    : int | None = None,
        phi_res      : int | None = None,    
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
        original_xs : np.ndarray
            array of original input Cartesian x coordinates.

        original_ys : np.ndarray
            array of original input Cartesian y coordinates.

        original_zs : np.ndarray
            array of original input Cartesian z coordinates.
        
        theta_range : tuple
            the polar (colatitude) range.
        
        phi_range : tuple
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

        # thetas are inclusive
        self._mercator_thetas = np.linspace(
            self._theta_range[0], 
            self._theta_range[1], 
            self._theta_res, dtype=_FLOAT
        )
        # phis are inclusive (only here if partial spherical grid)
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

        # generating a simple rectangular grid
        log.debug("Calling alloc_rect_map from C")
        self._map = platerecipy_clib_grid.alloc_rect_map(
            ctypes.c_int32(self._theta_res),
            ctypes.c_int32(self._phi_res),
            ctypes.c_int32(GTYPE_PARTIAL_SPHERICAL)
        )
        # assigning indices 
        npy_idxs = np.linspace(0, self._xs.size-1, self._xs.size, dtype=_INT)
        log.debug("Calling set_npy_idx from C")
        platerecipy_clib_grid.set_npy_idx(
            self._map,
            npy_idxs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        )

        # generating a pyvista mesh
        points = np.vstack([self._xs.ravel(), self._ys.ravel(), self._zs.ravel()]).T
        cells = np.zeros(5*(theta_res-1)*(phi_res-1), dtype=_INT)
        log.debug("Calling gen_cells_from_map from C")
        platerecipy_clib_grid.gen_cells_from_map(
            self._map, 
            cells.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        )
        cell_types = np.array([pv.CellType.QUAD for _ in range((theta_res-1)*(phi_res-1))])
        self._mesh = pv.UnstructuredGrid(cells, cell_types, points)

        log.debug('... grid created.')
    

    @property
    def theta_range(self) -> tuple[float, float]:
        '''Range of theta (polar angle)'''
        return self._theta_range

    @property
    def phi_range(self) -> tuple[float, float]:
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

    def __str__(self): # this will raise an error as interp_thete_res is not defined
        return f'A PartialSphericalGrid object with\ntheta_range={self.theta_range}\nphi_range={self.phi_range}\ninterp_theta_res={self.interp_theta_res}\ninterp_phi_res={self.interp_phi_res}\n'


    def interpolate_field(
        self,
        field           : np.ndarray,
        take_log        = False,
        method          = "tangent-plane", #"lat-lon", #"lat-lon" or "tangent-plane",
        **method_kwargs
    ) -> np.ndarray:
        """
        Interpolates the input field to the grid objects nodes.

        Parameters
        ----------
        field : np.ndarray
            array corresponding to the original Cartesian coordinate inputs.
        
        take_log : bool, False
            whether to perform the interpolation on the log space (suitable for
            fields that vary by orders of magnitude), default = False.
        
        method : str, default = "lat-lon"
            Whether to interpolate on the lat-lon space or to interpolate the 
            data points linearly on the tangent plane. Options: "lat-lon" or "tangent-plane"
        
        method_kwargs : dict, optional
            Additional arguments for `method = "tangent-plane"`. Namely, `k`, `eps`, 
            and `near_thresh`.
        
        Note
        ----
        `method = "tangent-plane"` is recommended for speed and accuracy.
        
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
        
        else:
            raise ValueError("Interpolation method is not recognized.")



        # forcing boundary conditions by averaging
        if np.abs(self.theta_range[0]) < 1e-3:
            # if close to the north pole
            grid_linear[ 0,  :] = np.mean(grid_linear[:,  0])
        if np.abs(self.theta_range[1] - np.pi) < 1e-3:
            # if close to the south pole
            grid_linear[-1,  :] = np.mean(grid_linear[:, -1])
        if (np.abs(-np.pi - self.phi_range[0]) < 1e-3) \
            and (np.abs(np.pi - self.phi_range[1]) < 1e-3):
            # if azimuthal endpoints are too close to each other
            # this should not happen
            grid_linear[ :,  0] = 0.5 * (grid_linear[:, 0] + grid_linear[:, -1])
            grid_linear[ :, -1] = grid_linear[:, 0]
            log.warning('Azimuthal end points are too close to each other. Please use a full spherical representation instead.')
        
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
        take_log    = False,
        method      = "tangent-plane", #"lat-lon" or "tangent-plane",
        **method_kwargs
    ) -> np.ndarray:
        """
        Maps the spherical field back to the original input format by interpolating 
        the nodes using the nearest method.

        Parameters
        ----------
        field : np.ndarray,
            array corresponding to the original Cartesian coordinate inputs.

        take_log : bool, False,
            whether to perform the interpolation on the log space (suitable for
            fields that vary by orders of magnitude), and makes a substantive
            difference if `method == 'linear'`, default = False.
        
        method : str, default = "lat-lon"
            Whether to interpolate on the lat-lon space or to interpolate the 
            data points linearly on the tangent plane. Options: "lat-lon" or "tangent-plane"
        
        method_kwargs : dict, optional
            Additional arguments for `method = "tangent-plane"`. Namely, `k`, `eps`, 
            and `near_thresh`.
        
        Note
        ----
        `method = "tangent-plane"` is recommended for speed and accuracy.
        
        Returns
        -------
        np.ndarray
        """
        log.debug('Interpolating fields back to the original input array')

        field = field.ravel()

        # checking if interpolation needs to occur on the log space
        if take_log:
            field = np.log(field)

        if method == "lat-lon":
            log.debug('... lat-lon linear interpolation')

            # by default the nearest method is used for the points
            # outside the convex hull and an integer field (e.g., plate ID)
            org_nearest = griddata(
                (self.thetas.ravel(), self.phis.ravel()),
                field.ravel(),  
                self.original_points_in_theta_phi, 
                method='nearest'
            )

            if field.dtype == _INT:
                org_order = org_nearest.astype(_INT)
            else:
                # a float field will be interpolated linearly within
                # the convex hull
                org_linear = griddata(
                    (self.thetas.ravel(), self.phis.ravel()),
                    field.ravel(),  
                    self.original_points_in_theta_phi, 
                    method='linear'
                )
                org_linear[np.isnan(org_linear)] = org_nearest[np.isnan(org_linear)]
                org_order = org_linear
                
        
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
            if not hasattr(self, '_kdtree'):
                log.debug('... constructing a similarity tree')
                self._kdtree = cKDTree(self._cart_mat)
            
            if (not hasattr(self, '_neigh_dists')) or (not hasattr(self, '_neighs')):
                log.debug('... querying the tree for near neighbors')
                self._neigh_dists, self._neighs = \
                    self._kdtree.query(self._original_cart_mat, k=method_kwargs['k'])
            
            neighs     = self._cart_mat[self._neighs]
            neigh_vals = field[self._neighs] 

            if field.dtype == _INT:
                org_order = neigh_vals[:, 0]
            else:
                # unit vectors on the surface
                #r_hats = self._original_cart_mat / self.r
                theta_hats = np.vstack([
                    np.cos(self.original_thetas.ravel())*np.cos(self.original_phis.ravel()),
                    np.cos(self.original_thetas.ravel())*np.sin(self.original_phis.ravel()),
                    -np.sin(self.original_thetas.ravel())
                ]).T
                phi_hats = np.vstack([
                    -np.sin(self.original_phis.ravel()),
                    np.cos(self.original_phis.ravel()),
                    np.zeros_like(self.original_thetas.ravel())
                ]).T

                # constructing a linear model on the tangent space
                log.debug('... constructing a linear LSQ model on the tangent space')
                ps = np.array([
                    neighs[:, i, :] - self._original_cart_mat[:, :] \
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
                org_order = (np.linalg.solve(GTG, GT@B)[:, 0]).ravel() # this may fail

                # if a point is too close to a gride point, simply use its value
                org_order[
                    self._neigh_dists[:, 0]/self.r < method_kwargs['near_thresh']
                ] = fs[0][
                    self._neigh_dists[:, 0]/self.r < method_kwargs['near_thresh']
                ]
                

        else:
            raise ValueError("Interpolation method is not recognized.")


        if take_log:
            org_order = np.exp(org_order)

        return org_order



class SphericalGrid(PartialSphericalGrid):
    """
    Uniform full spherical grid.
    """

    def __init__(
        self,
        original_xs  : np.ndarray,
        original_ys  : np.ndarray,
        original_zs  : np.ndarray,
        theta_res    : float | None = None,
        phi_res      : float | None = None,    
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
        Grid.__init__(self)

        theta_range = [0., np.pi]
        phi_range   = [-np.pi, np.pi]

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

        # thetas are [0, pi] inclusive
        self._mercator_thetas = np.linspace(
            self._theta_range[0], 
            self._theta_range[1], 
            self._theta_res, dtype=_FLOAT
        )
        # phis are [0, 2pi) excluding the end point
        self._mercator_phis = np.linspace(
            self._phi_range[0], 
            self._phi_range[1], 
            self._phi_res, dtype=_FLOAT,
            endpoint=False
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

        self._map = platerecipy_clib_grid.alloc_sph_map(
            ctypes.c_int32(self._theta_res),
            ctypes.c_int32(self._phi_res)
        )
        ni, nj = self._xs.shape
        # have to roll because in map, the north pole is the last index
        npy_idxs = np.roll(np.linspace(nj-1, nj*(ni-1), (ni-2)*nj+2, dtype=_INT), -1)
        
        platerecipy_clib_grid.set_npy_idx(
            self._map,
            npy_idxs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        )
        
        platerecipy_clib_grid.set_spherical_edge_length(
            self._map,
            self._xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._zs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_double(self._r)
        )

        # generating a pyvista mesh
        points = np.vstack(
            [self._xs.ravel(), self._ys.ravel(), self._zs.ravel()]
        ).T[phi_res-1:phi_res*(theta_res-1)+1, :]
        cells = np.zeros(
            5*(theta_res-3)*(phi_res) + 2*4*(phi_res), 
            dtype=_INT
        )
        log.debug("Calling gen_cells_from_map from C")
        platerecipy_clib_grid.gen_cells_from_map(self._map, cells.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        cell_types = np.array([
            pv.CellType.QUAD for _ in range((theta_res-3)*(phi_res))
        ] + [
            pv.CellType.TRIANGLE for _ in range((2*phi_res))
        ])
        self._mesh = pv.UnstructuredGrid(cells, cell_types, points)
        

    def enforce_data_consistency(self, array: np.ndarray):
        array[0, :] = array[0, -1]      # north pole
        #array[:, -1] = array[:, 0]      # wrap-around
        # wraparound column is removed ^^
        array[-1, :] = array[-1, 0]     # south-pole

    def __str__(self):
        return f'A SphericalGrid object with\ninterp_theta_res={self.interp_theta_res}\ninterp_phi_res={self.interp_phi_res}\n'
