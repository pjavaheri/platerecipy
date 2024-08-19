"""
@file model.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for model object and functions.
"""

import numpy as np
from scipy import ndimage
from skimage.segmentation import watershed

from inspect import isfunction

from platerecipy.transform import sph_fused_distance_threshold_transform as sphfdtt

def _tile_sph_surface(
    field           : np.ndarray,
    longitude_axis  : int
):
    if longitude_axis == 0:
        # the first index wraps around the sphere
        ni, nj = field.shape

        field_ext = np.zeros((3*ni, 3*nj), dtype=field.dtype)
        
        field_ext[     :   ni,      :   nj] = field[:, ::-1]
        field_ext[  ni : 2*ni,      :   nj] = field[:, ::-1]
        field_ext[2*ni : 3*ni,      :   nj] = field[:, ::-1]

        field_ext[     :   ni,   nj : 2*nj] = field[:, :]
        field_ext[  ni : 2*ni,   nj : 2*nj] = field[:, :]
        field_ext[2*ni : 3*ni,   nj : 2*nj] = field[:, :]

        field_ext[     :   ni, 2*nj : 3*nj] = field[:, ::-1]
        field_ext[  ni : 2*ni, 2*nj : 3*nj] = field[:, ::-1]
        field_ext[2*ni : 3*ni, 2*nj : 3*nj] = field[:, ::-1]

        field = field_ext

    elif longitude_axis == 1:
        # the second index wraps around the sphere
        ni, nj = field.shape

        field_ext = np.zeros((3*ni, 3*nj), dtype=field.dtype)
        
        field_ext[     :   ni,      :   nj] = field[::-1, :]
        field_ext[     :   ni,   nj : 2*nj] = field[::-1, :]
        field_ext[     :   ni, 2*nj : 3*nj] = field[::-1, :]
        
        field_ext[  ni : 2*ni,      :   nj] = field[:, :]
        field_ext[  ni : 2*ni,   nj : 2*nj] = field[:, :]
        field_ext[  ni : 2*ni, 2*nj : 3*nj] = field[:, :]
        
        field_ext[2*ni : 3*ni,      :   nj] = field[::-1, :]
        field_ext[2*ni : 3*ni,   nj : 2*nj] = field[::-1, :]
        field_ext[2*ni : 3*ni, 2*nj : 3*nj] = field[::-1, :]

        field = field_ext
        
    else:
        raise ValueError("unification_axis can be either None, 0, or 1.")

    return field_ext


class PlateModel(object):
    def __init__(self) -> None:
        self.stacked_field = None
        self.stack_weight_sum = 0.
        self._stacked_field_is_normalized = False

    def stack_field(
        self,
        field               : np.ndarray,
        invert              = False,
        take_log            = False,
        custom_function     = None,
        stack_weight        = 1.
    ) -> None:
        # a new field cannot be stacked if the stack is already normalized
        # (i.e., `find_plates` method has been called already)
        if self._stacked_field_is_normalized:
            raise Exception(
                "A new field cannot be stacked if the stack has been already \
                normalized. Ensure all necessary fields are stacked prior to \
                calling `find_plates`."
            )

        # adding the partial weight for normalization
        self.stack_weight_sum += stack_weight
        self._stacked_field_is_normalized = False

        # taking the log for fields that change by orders of magnitude
        if take_log:
            field = np.log10(field)
        
        # normalizing the field to conform to [0,1] range
        field = (field - field.min()) / (field.max() - field.min())

        # ensuring highs represent plate boundaries
        if invert:
            field = 1. - field
        
        # applying a custom transformation on the array
        if custom_function is not None:
            if isfunction(custom_function):
                field = custom_function(field)
            else:
                raise ValueError("The custom_function must be a proper function.")

        # if stack is empty
        if self.stacked_field is None:
            self.stacked_field = field
        else:
            self.stacked_field += field * stack_weight
    
    def _normalize_stacked_field(self):
        if not self._stacked_field_is_normalized:
            self.stacked_field *= (1./self.stack_weight_sum)
            self._stacked_field_is_normalized = True

    
    def find_plates(
        self,
        longitudinal_axis       = None,
        interior_quantile       = 0.,
        boundary_quantile       = 0.9,
        spatial_tolerance       = None,
        spatial_weight          = None,
        num_threads             = 1,
        min_marker_size         = None,
        watershed_connectivity  = 1,
        watershed_compactness   = 1.,
        halo_quantile           = 1.,
        halo_spatial_tolerance  = None
    ) -> np.ndarray:
        """
        Applies segmentation on `stacked_field` and returns an integer array of 
        the same shape with each cell carrying the plate ID (i.e., segment
        number). This function can be called multiple times from the same object
        instance as it does not alter 
        
        Warning
        -------
        There are several options that fundamentally alter the way this function 
        works. The user is suggested to carefully read the documentation for all 
        input parameters regardless of their apparent relevance to the specific 
        task at hand.

        Parameters
        ----------
        longitudinal_axis : int, optional
            Which axis (0 or 1) represents the longitude for the mercator 
            projection. By inputting a number, the operations are altered
            knowing `stacked_field` is a mercator projection of a spherical
            surface with `longitudinal_axis` representing the axis along which
            the field wraps around the sphere (the other axis points toward the 
            poles). 

            In consequence, instead of a normal Euclidean distance transform, 
            a spherical distance transform is used. Additionally, plate IDs will
            we adjusted to ensure boundary conditions for a mercator projection 
            (e.g., the plates that wrap around have the same ID).

        interior_quantile : float, default=0.
            The quantile that is sure to be a plate interior. Consequently, all
            values less than the value that corresponds to this quantile will be
            replaced with 0. It is suggested to leave this option at 0 (default).
        
        boundary_quantile : float, default=0.9
            The quantile that represents values above which that are presumed to 

        spatial_tolerance : float, optional
            Whether to use a distance transform and consider regions with 
            distances less than `spatial_tolerance` from the boundary (determined
            by `boundary_quantile`) as a part of the boundary as well. This 
            option is useful when segments (plates) have imperfect boundaries 
            that require some spatial tolerance to close them.

            If `longitudinal_axis` is provided, it is presumed that the field 
            is a mercator projection and instead of a Euclidean distance 
            transform, a spherical one will be applied. In that case, 
            `spatial_tolerance` will be considered as the angle of tolerance 
            on the great circle passing through a given pair of points in radians.

        spatial_weight          = None,
        num_threads             = 1,
        min_marker_size         = None,
        watershed_connectivity  = 1,
        watershed_compactness   = 1.,
        halo_quantile           = 1.,
        halo_spatial_tolerance  = None

        """
        # normalizing the stacked fields
        self._normalize_stacked_field()
        
        
        self.boundary_quantile_value, \
            self.interior_quantile_value, \
            self.halo_quantile_value = np.quantile(
            self.stacked_field, 
            [boundary_quantile, interior_quantile, halo_quantile]
        )

        # operating on a separate copy
        self.stacked_field_for_watershed = self.stacked_field.copy()

        # forcing specific regions to be an interior, and thus, show up as a marker
        # it is suggested not to be used!
        if interior_quantile > 0.:
            self.stacked_field_for_watershed[
                self.stacked_field_for_watershed < self.interior_quantile_value
            ] = 0.

        self.markers = self.stacked_field_for_watershed < self.boundary_quantile_value
        
        if longitudinal_axis is not None:
            # field is spherical
            if not ((longitudinal_axis == 0) or (longitudinal_axis == 1)):
                raise ValueError("Longitudinal axis must be either 0 or 1.")
            
            if halo_spatial_tolerance is not None:
                halo_region = ~sphfdtt(
                    self.stacked_field_for_watershed > self.halo_quantile_value, 
                    halo_spatial_tolerance,
                    num_threads=num_threads
                )
                self.stacked_field_for_watershed[halo_region] = self.halo_quantile_value

            if spatial_tolerance is not None:
                complement_markers = ~sphfdtt(
                    self.stacked_field_for_watershed > self.boundary_quantile_value, 
                    spatial_tolerance,
                    num_threads=num_threads
                )
        else:
            # field is Cartesian
            if halo_spatial_tolerance is not None:
                halo_region = ndimage.distance_transform_edt(
                    self.stacked_field_for_watershed < self.halo_quantile_value
                ) > halo_spatial_tolerance
            self.stacked_field_for_watershed[halo_region] = self.halo_quantile_value

            if spatial_tolerance is not None:
                complement_markers = ndimage.distance_transform_edt(
                    self.stacked_field_for_watershed < self.boundary_quantile_value
                ) > spatial_tolerance
        
        if spatial_tolerance is not None:

            # whether to also elevate the region near the boundary
            # it is suggested not to be used! 
            if spatial_weight is not None:
                self.stacked_field_for_watershed += spatial_weight*(
                        (~complement_markers).astype(self.stacked_field.dtype)
                    )
                self.stacked_field_for_watershed *= 1./(1. + spatial_weight)

            temp_labels = ndimage.label(self.markers)[0]

            # removing labels that have non-empty intersections
            temp_unique_labels = np.unique(temp_labels[temp_labels != 0])

            for label in temp_unique_labels:
                if np.all(complement_markers[temp_labels == label] == False):
                    # a new valid region is found which was overwhelmed by the
                    # spatial tolerance
                    complement_markers[temp_labels == label] = True
            
            del temp_labels, temp_unique_labels
            self.markers = complement_markers

        # filtering out micro markers
        if min_marker_size is not None:
            temp_labels = ndimage.label(self.markers)[0]

            # removing labels that have non-empty intersections
            temp_unique_labels = np.unique(temp_labels)

            for label in temp_unique_labels:
                if np.sum(temp_labels == label) < min_marker_size:
                    self.markers[temp_labels == label] = False
            
            del temp_labels, temp_unique_labels
        
        # converting to an 8-bit integer (i.e., pixel-like) field
        field_uint8 = (255.*self.stacked_field_for_watershed).astype(np.uint8)

        if longitudinal_axis is not None:
            # using tiling operation to ensure spherical boundary conditions
            markers_ext = _tile_sph_surface(
                field=self.markers, 
                longitude_axis=longitudinal_axis
            )

            field_uint8_ext = _tile_sph_surface(
                field=field_uint8, 
                longitude_axis=longitudinal_axis
            )

            labels_ext = ndimage.label(markers_ext)[0]
            plates_ext = watershed(
                field_uint8_ext, 
                labels_ext, 
                connectivity=watershed_compactness, 
                compactness=watershed_connectivity
            )

            ni, nj = self.stacked_field.shape
            self.plate_IDs = plates_ext[ni:2*ni, nj:2*nj]

            if longitudinal_axis == 0:
                # fixing the polar condition
                pole1_IDs = np.unique(self.plate_IDs[:, 0])
                for ID in pole1_IDs:
                    self.plate_IDs[self.plate_IDs == ID] = pole1_IDs.min()
                pole2_IDs = np.unique(self.plate_IDs[:, -1])
                for ID in pole2_IDs:
                    self.plate_IDs[self.plate_IDs == ID] = pole2_IDs.min()

                # fixing the wrapping boundaries
                side1_IDs = self.plate_IDs[0, :]
                side2_IDs = self.plate_IDs[-1, :]

                for i in range(self.plate_IDs.shape[1]):
                    if side1_IDs[i] != side2_IDs[i]:
                        self.plate_IDs[self.plate_IDs == side2_IDs[i]] = side1_IDs[i]

                raw_unified_IDs = np.unique(self.plate_IDs)
                for i, ID in enumerate(raw_unified_IDs):
                    self.plate_IDs[self.plate_IDs == ID] = -i
                self.plate_IDs *= -1
                
            else:
                # fixing the polar condition
                pole1_IDs = np.unique(self.plate_IDs[0, :])
                for ID in pole1_IDs:
                    self.plate_IDs[self.plate_IDs == ID] = pole1_IDs.min()
                pole2_IDs = np.unique(self.plate_IDs[-1, :])
                for ID in pole2_IDs:
                    self.plate_IDs[self.plate_IDs == ID] = pole2_IDs.min()

                # fixing the wrapping boundaries
                side1_IDs = self.plate_IDs[:, 0]
                side2_IDs = self.plate_IDs[:, -1]

                for i in range(self.plate_IDs.shape[0]):
                    if side1_IDs[i] != side2_IDs[i]:
                        self.plate_IDs[self.plate_IDs == side2_IDs[i]] = side1_IDs[i]

                raw_unified_IDs = np.unique(self.plate_IDs)
                for i, ID in enumerate(raw_unified_IDs):
                    self.plate_IDs[self.plate_IDs == ID] = -i
                self.plate_IDs *= -1
        else:
            # it is a Cartesian plane
            labels = ndimage.label(self.markers)[0]
            self.plate_IDs = watershed(
                field_uint8, 
                labels,
                connectivity=watershed_compactness, 
                compactness=watershed_connectivity
            )
            
        return self.plate_IDs.copy()
