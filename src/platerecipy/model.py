"""
File brief
----------
`model.py`

Module for model object and functions.

This is a part of `platerecipy` package. For license and citation, please
refer to the main repository:
[github.com/pjavaheri/platerecipy](github.com/pjavaheri/platerecipy)

Author(s): 
Pejvak Javaheri; [pejvak.javaheri@mail.utoronto.ca](mailto:pejvak.javaheri@mail.utoronto.ca)
"""

import logging
log = logging.getLogger(__name__)

import numpy as np
from scipy.spatial import cKDTree

from .transform import label_markers, distance_transform
from .segmentation import random_walker
from .grid import Grid
from . import _FLOAT


def _dict_arg_max(d):
    keys = list(d.keys())
    maxval = d[keys[0]]
    maxkey = int(keys[0])
    for key in keys:
        if d[key] > maxval:
            maxkey = int(key)
    return maxkey

def _get_pointwise_omega(grid: Grid, velocity_key: str):
    omega_hat = np.cross(grid.mesh.points, grid.mesh[velocity_key])
    omega_norm = np.linalg.norm(omega_hat, axis=1)
    omega_hat /= omega_norm[:, None]
    omega_mag = np.linalg.norm(grid.mesh[velocity_key], axis=1)/grid.r
    return omega_mag, omega_hat

def _evolve_model(grid: Grid, omega_mag: np.ndarray, omega_hat: np.ndarray, dt: float):
    angle = omega_mag*dt
    new_pos = grid.mesh.points * (np.cos(angle)[:, None]) \
            + np.cross(omega_hat, grid.mesh.points)*(np.sin(angle)[:, None]) \
            + omega_hat * ((np.vecdot(omega_hat, grid.mesh.points)*(1-np.cos(angle)))[:, None])
    return new_pos


class PlateModel(object):
    """
    A class for defining plate detection parameters, field stacking, 
    and performing segmentation.
    """

    def __init__(self, grid: Grid) -> None:
        """
        Create a `PlateModel` object by initializing it using a `Grid`.

        Parameters
        ----------
        grid : Grid,
            An instance of a grid 
        """
        self.stacked_field = None
        self.stack_weight_sum = 0.
        self._stacked_field_is_normalized = False

        self.grid = grid

        self.frwd_prop_IDs = None
        self.bkwd_prop_IDs = None

    
    def clear_stacked_field(self) -> None:
        """
        To reset and clear the stacked field.
        """
        self.stacked_field = None
        self.stack_weight_sum = 0.
        self._stacked_field_is_normalized = False
        

    def stack_field(
        self,
        field               : np.ndarray,
        invert              = False,
        take_log            = False,
        stack_weight        = 1.
    ) -> None:
        """
        Stack a new field by normalizing the field values.

        Parameters
        ----------
        field : np.ndarray,
            The field to be stacked.
            
        invert : bool, default=False,
            If increase in input `field` corresponds to a decrease in deformation.
        
        take_log : bool, default=False,
            If the input `field` varies by orders of magnitude.

        stack_weight : float, default=1.,
            The corresponding weight when stacked (between 0 and 1).
        """

        # a new field cannot be stacked if the stack is already normalized
        # (i.e., `find_plates` method has been called already)
        if self._stacked_field_is_normalized:
            raise Exception(
                "A new field cannot be stacked if the stack has been already \
                normalized. Ensure all necessary fields are stacked prior to \
                calling `find_plates`."
            )
        
        # enforcing row-major structure
        field = field.astype(order='C', dtype=_FLOAT, copy=False)

        # adding the partial weight for normalization
        self.stack_weight_sum += stack_weight
        self._stacked_field_is_normalized = False

        # taking the log for fields that change by orders of magnitude
        if take_log:
            field = np.log(field)
        
        # normalizing the field to conform to [0,1] range
        field = (field - field.min()) / (field.max() - field.min())

        # ensuring highs represent plate boundaries
        if invert:
            field = 1. - field
                
        # if stack is empty
        if self.stacked_field is None:
            self.stacked_field = field * stack_weight
        else:
            self.stacked_field += field * stack_weight
    
    def _normalize_stacked_field(self):
        if not self._stacked_field_is_normalized:
            self.stacked_field *= (1./self.stack_weight_sum)
            self._stacked_field_is_normalized = True

    # a reset stacked field to be added
    
    def find_plates(
        self,
        boundary_quantile       : float = 0.9,
        boundary_absolute       : float = 1.0,
        separation_tolerance    : float | None = None,
        num_threads             : int = 1,
        min_marker_size         : int | None = None,
        preserve_small_markers  : bool = False,
        manual_markers          : np.ndarray | None = None,
        identify_nonconforming  : bool = False,
        RW_beta                 : float = 100.,
        RW_solver_tolerance     : float = 1e-3,
        RW_solver               : str = 'LU',
        return_IDs              : bool = True
    ) -> np.ndarray | None:
        """
        Applies segmentation on `stacked_field` and returns an integer array of 
        the same shape with each cell carrying the plate ID (i.e., segment
        number). This function can be called multiple times from the same object
        instance as it does not alter `stacked_field`.
        
        Warning
        -------
        There are several options that fundamentally alter the way this function 
        works. The user is suggested to carefully read the documentation for all 
        input parameters regardless of their apparent relevance to the specific 
        task at hand.

        Warning
        -------
        If both `boundary_quantile` and `boundary_absolute` are specified, a 
        node qualifies as a possible boundary point so long as it satisfies one 
        of the two conditions.

        Parameters
        ----------
        boundary_quantile : float, default=0.9
            The quantile that represents values above which that are presumed to 
            feature plate boundaries.
        
        boundary_absolute : float, default=1.0
            The absolute value that represents values above which that are 
            presumed to feature plate boundaries.

        separation_tolerance : float, optional
            Whether to use a distance transform and consider regions with 
            distances less than `separation_tolerance` from the boundary (determined
            by `boundary_quantile`) as a part of the boundary as well. This 
            option is useful when segments (plates) have imperfect boundaries 
            that require some separation tolerance to close them.

            For a `SphericalGrid`, it is presumed that the field 
            is a mercator projection and instead of a Euclidean distance 
            transform, a spherical one will be applied. In that case, 
            `separation_tolerance` will be considered as the angle of tolerance 
            on the great circle passing through a given pair of points in radians.

        num_threads : int, default=1
            Number of threads to use to perform the spherical distance transform.
        
        min_marker_size : int, optional
            If provided, watershed markers will be filtered such that markers 
            with fewer cells/pixels than `min_marker_size` will be ignored. This
            is useful when then input field is noisy or not coherent enough.

        preserve_small_markers : bool, default=False
            Whether to reinstate small markers obscured by by the separation
            tolerance step.
        
        manual_markers : np.ndarray, optional
            This option ignores all previous arguments relating to marker
            specification. Instead, it allows the user to manually provide 
            a 2D Boolean array of the same shape as `self.stacked_field` that 
            indicates the seeds to the RW algorithm.

        identify_non_conforming : bool, default=False
            Whether to extract and separately label non-conforming regions 
            defined by featuring a stacked field greater than 0.5 and a ID 
            probability less than 0.5 with ID=0.

        RW_beta : float, default=100.
            Gaussian beta parameter for random walker connection weights. This 
            parameter controls the sharpness of the boundaries.
        
        RW_solver_tolerance : float, default=1e-3
            Tolerance value for the choice of solver that requires a tolerance
            (e.g., if `RW_solver = 'CG'`).
        
        RW_solver : str, default='LU'
            Numerical solver used to obtain an random walker probability solution.
            Possible choices: 'direct', 'LU', 'CG', or 'FA'.
        
        return_IDs : bool, default=True
            Whether to return a copy of `self.plate_IDs`.
        """
        if self.stacked_field is None:
            raise ValueError('No field is stacked for segmentation. Please use `stack_field()` method function to stack at least one field.')

        # normalizing the stacked fields
        self._normalize_stacked_field()
        
        if manual_markers is not None:
            self.markers = manual_markers
        else:

            self.boundary_quantile_value = np.quantile(
                self.stacked_field, 
                [boundary_quantile]
            )[0]

            self.boundary_absolute = boundary_absolute

            self.boundary_threshold = np.min([self.boundary_absolute, self.boundary_quantile_value])

            # operating on a separate copy
            #self.stacked_field_for_segmentation = self.stacked_field.copy()
            
            self.markers = (self.stacked_field < self.boundary_threshold)
            
            # filtering out micro markers
            if min_marker_size is not None:
                log.info(f"Filtering marker pataches comprised of fewer than {min_marker_size} nodes ...")
                """
                temp_labels = ndimage.label(self.markers)[0]
                """
                temp_labels = label_markers(self.grid, self.markers)

                # removing labels that have non-empty intersections
                temp_unique_labels = np.unique(temp_labels)

                for label in temp_unique_labels:
                    if np.sum(temp_labels == label) < min_marker_size:
                        self.markers[temp_labels == label] = False
                
                del temp_labels, temp_unique_labels
            
            if separation_tolerance is not None:
                log.info(f"Allowing for {separation_tolerance:5.2e} (in radians) of separation tolerance to patch missing plate boundary segments ...")
                if preserve_small_markers:
                    log.info(f"Preserving small markers against separation tolerance ...")
                    temp_pre_labels = label_markers(self.grid, self.markers)
            
                self.old_markers = self.markers.copy()
                self.markers = distance_transform(
                    grid        = self.grid, 
                    markers     = self.markers, 
                    threshold   = separation_tolerance, 
                    num_threads = num_threads
                )
            
                if preserve_small_markers:
                    temp_unique_labels = np.unique(temp_pre_labels)
                    temp_unique_labels = temp_unique_labels[temp_unique_labels>0]
                    
                    for label in temp_unique_labels:
                        temp_mask = (temp_pre_labels == label)
                        if not np.any(self.markers[temp_mask]):
                            self.markers[temp_mask] = True

        # labeling the markers with positive integers
        self.labels = label_markers(self.grid, self.markers)

        # handling both spherical and planar cases
        log.info(f"Applying the Random Walker algorithm with beta={RW_beta:5.2e} ...")
        self.plate_IDs, self.ID_probs = random_walker(
            data             = self.stacked_field, #self.stacked_field_for_segmentation,
            labels           = self.labels,
            beta             = RW_beta,
            solver_tol       = RW_solver_tolerance,
            solver           = RW_solver,
            grid             = self.grid
        )
        
        if identify_nonconforming:
            self.plate_IDs[(self.stacked_field > self.boundary_threshold) | (self.ID_probs < 0.5)] = 0
        
        if return_IDs:
            return self.plate_IDs.copy()