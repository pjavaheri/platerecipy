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
from scipy import ndimage

from .transform import gridded_fused_distance_threshold_transform, \
                        pgridded_fused_distance_threshold_transform
from .segmentation import random_walker
from .grid import Grid, PartialSphericalGrid, SphericalGrid
from . import _FLOAT


def _unify_wraparound_labels(labels: np.ndarray) -> np.ndarray:
    """
    Ensures marker labels conform to the azimuthal continuity (i.e., wraparound
    boundary condition).

    Parameters
    ----------
    labels : np.ndarray,
        A 2D integer array with labels positive integers and unmarked regions 
        zeros.

    Returns
    -------
    np.ndarray

    Warning
    -------
    The input array, `labels`, will be modified, but the modification is not 
    guaranteed to be sequential. `_make_labels_sequential()` should be called 
    subsequently.

    Warning
    -------
    It is assumed that the first dimension corresponds to theta (polar angle) 
    and the second to phi (azimuthal angle).
    """
    for i in range(labels.shape[0]):
        if labels[i, 0] > 0 and labels[i, -1] > 0:
            labels[labels == labels[i, -1]] = labels[i, 0]
    
    return labels

def _make_labels_sequential(labels: np.ndarray) -> np.ndarray:
    """
    Ensures marker labels are sequential and no gaps between positive IDs.

    Parameters
    ----------
    labels : np.ndarray,
        A 2D integer array with labels positive integers and unmarked regions 
        zeros.

    Returns
    -------
    np.ndarray

    Warning
    -------
    The input array, `labels`, will be modified, but the modification is not 
    guaranteed to be sequential. `_make_labels_sequential()` should be called 
    subsequently.
    
    Warning
    -------
    It is assumed that the first dimension corresponds to theta (polar angle) 
    and the second to phi (azimuthal angle).
    """
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels != 0]
    
    for i in range(unique_labels.size):
        labels[labels==unique_labels[i]] = i+1
    
    return labels



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
        # if a spherical grid, wraparound azimuthally 
        self._wraparound_azimuthally = isinstance(grid, SphericalGrid)
        self._use_spherical_distance = isinstance(grid, PartialSphericalGrid)
    
    @property
    def wraparound_azimuthally(self):
        '''Whether to apply wraparound boundary conditions along phi (i.e., for a SphericalGrid).'''
        return self._wraparound_azimuthally

    @property
    def use_spherical_distance(self):
        '''Whether to apply use the great circle angle of separation instead of planar distance transform (i.e., for a SphericalGrid).'''
        return self._use_spherical_distance
    
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
        boundary_quantile       = 0.9,
        boundary_absolute       = 1.0,
        separation_tolerance    = None,
        num_threads             = 1,
        min_marker_size         = None,
        preserve_small_markers  = False,
        manual_markers          = None,
        identify_nonconforming  = False,
        RW_beta                 = 100.,
        RW_solver_tolerance     = 1e-3,
        RW_solver               = 'LU',
        return_IDs              = True
    ) -> np.ndarray:
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
            )

            self.boundary_absolute = boundary_absolute

            # operating on a separate copy
            self.stacked_field_for_segmentation = self.stacked_field.copy()
            
            self.markers = (self.stacked_field_for_segmentation < self.boundary_quantile_value) \
                & (self.stacked_field_for_segmentation < self.boundary_absolute)
            
            # filtering out micro markers
            if min_marker_size is not None:
                temp_labels = ndimage.label(self.markers)[0]

                # removing labels that have non-empty intersections
                temp_unique_labels = np.unique(temp_labels)

                for label in temp_unique_labels:
                    if np.sum(temp_labels == label) < min_marker_size:
                        self.markers[temp_labels == label] = False
                
                del temp_labels, temp_unique_labels
            
            if preserve_small_markers:
                temp_pre_labels = ndimage.label(self.markers)[0]
                if self.wraparound_azimuthally:
                    _unify_wraparound_labels(temp_pre_labels)
            
            if separation_tolerance is not None:
                if self.use_spherical_distance and self.wraparound_azimuthally:
                    # full sphere
                    # separation_tolerance is treated as radians on the great-circle
                    self.markers = ~gridded_fused_distance_threshold_transform(
                        xs          = self.grid.xs, 
                        ys          = self.grid.ys, 
                        zs          = self.grid.zs,
                        arr         = ~self.markers, 
                        R           = self.grid.r,
                        threshold   = separation_tolerance,
                        num_threads = num_threads
                    )
                elif self.use_spherical_distance:
                    # a partial sphere
                    self.markers = ~pgridded_fused_distance_threshold_transform(
                        xs          = self.grid.xs, 
                        ys          = self.grid.ys, 
                        zs          = self.grid.zs,
                        theta_range = self.grid.theta_range,
                        phi_range   = self.grid.phi_range,
                        arr         = ~self.markers, 
                        R           = self.grid.r,
                        threshold   = separation_tolerance,
                        num_threads = num_threads
                    )
                else: 
                    # a flat grid
                    pass
                    # planar distance transform with the separation_tolerance treated 
                    # as the distance in terms of unit grid spacing
                    self.markers = ~(
                        ndimage.distance_transform_edt(
                            input=~self.markers
                        ) < separation_tolerance
                    )
            
            if preserve_small_markers:
                temp_unique_labels = np.unique(temp_pre_labels)
                temp_unique_labels = temp_unique_labels[temp_unique_labels>0]
                
                for label in temp_unique_labels:
                    temp_mask = (temp_pre_labels == label)
                    if not np.any(self.markers[temp_mask]):
                        self.markers[temp_mask] = True

        # labeling the markers with positive integers
        labels, _ = ndimage.label(self.markers)
        
        if self.wraparound_azimuthally:
            # making labels consistent and sequential
            _unify_wraparound_labels(labels)
            _make_labels_sequential(labels)
        
        # handling both spherical and planar cases

        self.plate_IDs, self.ID_probs = random_walker(
            data             = self.stacked_field_for_segmentation,
            labels           = labels,
            beta             = RW_beta,
            solver_tol       = RW_solver_tolerance,
            solver           = RW_solver,
            grid             = self.grid
        )
        
        if identify_nonconforming:
            self.plate_IDs[(self.stacked_field > 0.5) | (self.ID_probs < 0.5)] = 0
        
        if return_IDs:
            return self.plate_IDs.copy()