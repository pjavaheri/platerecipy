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

    def stack_field(
        self,
        field               : np.ndarray,
        invert              = False,
        take_log            = False,
        custom_function     = None,
        stack_weight        = 1.
    ) -> None:
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
    
    def find_plates(
        self,
        longitudinal_axis       = None,
        interior_quantile       = None,
        boundary_quantile       = None,
        spatial_tolerance       = 0.,
        spatial_weight          = 0.,
        num_threads             = 1,
        min_marker_size         = None,
        watershed_connectivity  = 1,
        watershed_compactness   = 1.
    ) -> np.ndarray:
        
        self.boundary_quantile_value, self.interior_quantile_value = np.quantile(
            self.stacked_field, 
            [boundary_quantile, interior_quantile]
        )

        self.stacked_field_for_watershed = self.stacked_field.copy()

        #self.stacked_field_for_watershed[
        #    self.stacked_field_for_watershed < self.interior_quantile_value
        #] = 0.

        self.markers = self.stacked_field_for_watershed < self.boundary_quantile_value
        
        if longitudinal_axis is not None:
            # field is spherical
            if not ((longitudinal_axis == 0) or (longitudinal_axis == 1)):
                raise ValueError("Longitudinal axis must be either 0 or 1.")
            if spatial_tolerance > 0.:
                complement_markers = ~sphfdtt(
                    self.stacked_field_for_watershed > self.boundary_quantile_value, 
                    spatial_tolerance,
                    num_threads=num_threads
                )
        else:
            # field is Cartesian
            if spatial_tolerance > 0.:
                complement_markers = ndimage.distance_transform_edt(
                    self.stacked_field_for_watershed < self.boundary_quantile_value
                ) > spatial_tolerance
        
        if spatial_tolerance > 0.:
            self.stacked_field_for_watershed = self.stacked_field_for_watershed \
                + spatial_weight*(
                    (~complement_markers).astype(self.stacked_field.dtype)
                )
            self.stacked_field_for_watershed /= (1. + spatial_weight)

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
            
        self.markers[
            self.stacked_field < self.interior_quantile_value
        ] = True

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
        field_uint8 = np.uint8(255.*self.stacked_field_for_watershed)

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
