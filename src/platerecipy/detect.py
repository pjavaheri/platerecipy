"""
@file detect.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for detection functions.
"""

import numpy as np

from scipy import ndimage as ndi

from skimage.morphology import disk
from skimage.segmentation import watershed
from skimage.filters import rank
from skimage.util import img_as_ubyte

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

def get_plate_IDs(
    field               : np.ndarray,
    invert              = False,
    take_log            = False,
    spatial_tolerance   = 0.,
    spatial_weight      = 1.,
    boundary_quantile   = 0.95,
    interior_quantile   = 0.5,
    feature_thickness   = None,
    max_num_plates      = None,
    min_plate_size      = None,
    kernel              = 'gradient',
    kernel_threshold    = None,
    marker_threshold    = None,
    unification_axis    = None,
    interior_cutoff     = None,
    diagnostics         = False
) -> np.ndarray:
    """
    Identifies plates and returns an integer field representing each cell's 
    plate ID.

    Warning
    -------
    The values contained in `field` must be such that either the highs or the 
    lows represent plate boundaries. Meaning, a field such as raw surface 
    temperature is not suitable as median values tend to mostly comprise plate 
    interiors (instead of one extreme).
    In such cases, for example, `field` can be:
        >>> field = absolute_value(bad_field - mean(bad_field))
        >>> IDs = get_plate_IDs(field, ...)


    Parameters
    ----------
    field: np.ndarray
        A 2D scalar array containing values that are used to determine plate 
        interior and boundaries. There are no restrictions on the array type 
        (float | int | bool), values (positive | negative), or order (linear | 
        exponential). 
        However, other parameters must be set/changed accordingly to ensure a 
        sensible results.
    
    invert: bool, optional
        By default, plate boundaries are considered to be the highs of the 
        field. Otherwise, `field` must be inverted. Meaning, only set 
        invert=True if the lows correspond to plate boundaries.
        Default: False  
        
    take_log: bool, optional
        Whether the values change by orders of magnitude such that plate 
        recognition is better performed on the logarithm of `field`, instead.
        Default: False
    
    feature_thickness: int, optional
        Approximate average thickness of plate features (e.g., plate boundary 
        thickness) in pixels. If not provided, the average thickness will be 
        automatically calculated by applying the distance transform and 
        calculating the average of non-zero values.
        Default: None
    
    max_num_plates: int, optional
        Maximum number of plates to be found. By default, there is no cap on the 
        number of plates. If set, markers will be rank ordered according to 
        their size and plates exceeding `max_num_plates` which are smaller in 
        marker surface area are rejected. It should be noted that the rank
        ordering is not exact as it is applied on the input to the watershed 
        function.
        Note that this will be applied after `min_plate_size`.
        Default: None
    
    min_plate_size: int, optional
        The minimum number of marker pixels for a region to become a candidate
        plate. Similar to `max_num_plates`, this is applied on the markers that
        are used to feed the watershed mechanism and do not perfectly map to 
        plate size.
        Note that this will be applied before `max_num_plates`.
        Default: None
    
    kernel: str, optional
        The mechanism to detect plate boundaries, and subsequently, plate 
        interiors to be used as basins to feed the watershed mechanism.

        The following kernels are available:
        * 'direct':
            The original field will be directly compared against 
            `kernel_threshold`.
            
        * 'gradient': 
            The gradient will be applied to `field` and plate boundaries will
            show up as surges in the gradient compared to the approximately 
            zero-gradient plate interiors. However, diffuse plate boundaries 
            that have thicknesses much larger than `feature_thickness` may be
            considered a separate plate since within a large diffuse boundary,
            gradient drops.
        
        * function:
            A function accepting one argument with type np.ndarray and similar
            shape as the input field. The output will be set compared against
            `kernel_threshold`.
        
        Default: 'gradient'

    kernel_threshold:


    marker_threshold: 

    unification_axis: int, optional
        Whether to enforce boundary conditions on plate IDs by unifying plate
        IDs along the periodic axis `unification_axis`. If specified, `field` 
        will be extended by repeating along `unification_axis` and mirroring 
        along the other. The value should be either 0 or 1 and should correspond 
        to the azimuthal (i.e., the longitudinal) direction that wraps around 
        the sphere. 
        Note that enabling this option brings additional benefits to plate
        recognition: approximity to the array boundaries (left | right | top |
        bottom) can bring errors to plate detection. The explicit tiling 
        eliminates (or significantly reduces) that issue.

    interior_cutoff: float, optional
        proportion (quantile) of the surface that is sure to be plate interior. 
        This value is necessary for when `kernel= 'distance'`. If not provided 
        and only if `kernel= 'distance'`, 0.25 will be assumed.
        Default: None (or 0.25 if `kernel= 'distance'`)

    diagnostics: bool, optional
        Whether to produce detailed step by step report of the procedures 
        applied.
    """
    # taking the log for fields that change by orders of magnitude
    if take_log:
        field = np.log10(field)
    
    # normalizing the field to conform to [0,1] range
    field = (field - field.min()) / (field.max() - field.min())

    # ensuring highs represent plate boundaries
    if invert:
        field = 1. - field
    
    boundary_quantile_value = np.quantile(field, boundary_quantile)
    interior_quantile_value = np.quantile(field, interior_quantile)

    if spatial_tolerance > 0.:
        field = (
            field
            + spatial_weight * sphfdtt(
                    field > boundary_quantile_value, 
                    spatial_tolerance
                ).astype(field.dtype)
        ) / (1. + spatial_weight)

    # converting to an 8-bit integer (i.e., pixel-like) field
    field = np.uint8(255*field)

    return field #test

    if feature_thickness is None:
        # thickness needs to be estimated
        feature_thickness = np.minimum(field.shape[0], field.shape[1]) // 16
        
    if unification_axis is not None:
        field_ext = _tile_sph_surface(field=field, longitude_axis=unification_axis)

    if kernel == 'direct':
        markers = field < kernel_threshold
    elif kernel == 'gradient':
        # gradient with a footprint
        markers = rank.gradient(
            field, 
            footprint=disk(feature_thickness // 2)
        ) < kernel_threshold
    
    elif kernel == 'laplacian':
        # bilateral sum of the Laplacian
        markers = rank.sum_bilateral(
            ndi.laplace(field),
            footprint=disk(feature_thickness // 2),
            s0=10,
            s1=10
        ) < kernel_threshold  # should this be 1?
        
        #elif isfunction(kernel):
        #markers = kernel(
        #    field
        #) < kernel_threshold
    
    elif isfunction(kernel):
        markers = kernel(
            field
        )
        

    else:
        # unrecognized kernel
        raise ValueError("Unrecognized kernel.")

    basin_labels = ndi.label(markers)[0]
    unique_labels = np.unique(basin_labels)
    print(f"Found {len(unique_labels)} raw markers for processing.")
    
    # filtering basin labels based on their size
    if (min_plate_size is not None) or (max_num_plates is not None):
        label_to_size = {}
        for label in unique_labels:
            if label != 0:
                label_to_size[label] = np.sum(basin_labels == label)
        print(label_to_size)
    
    if min_plate_size is not None:
        unique_labels_copy = unique_labels.copy()
        for label in unique_labels_copy:
            if (label != 0) and (label_to_size[label] < min_plate_size):
                unique_labels[unique_labels == label] = 0
                basin_labels[basin_labels == label] = 0
                label_to_size[label] = None
        unique_labels = np.unique(unique_labels)

    if max_num_plates is not None:
        if unification_axis is not None:
            max_num_plates *= 9
        labels_sorted = list(label_to_size.keys())
        labels_sorted.sort(reverse=True)
        
        for label in labels_sorted[max_num_plates:]:
            unique_labels[unique_labels == label] = 0
            basin_labels[basin_labels == label] = 0
            label_to_size[label] = None

    #if kernel == 'gradient':
        #finer_gradient = rank.gradient(field, disk(feature_thickness // 4))
        #plate_IDs = watershed(field, basin_labels)

    
    plate_IDs = watershed(field, basin_labels)

    if unification_axis is not None:
        plate_IDs = plate_IDs[ni:2*ni, nj:2*nj]

        if unification_axis == 0:
            # fixing the polar condition
            pole1_IDs = np.unique(plate_IDs[:, 0])
            for ID in pole1_IDs:
                plate_IDs[plate_IDs == ID] = pole1_IDs.min()
            pole2_IDs = np.unique(plate_IDs[:, -1])
            for ID in pole2_IDs:
                plate_IDs[plate_IDs == ID] = pole2_IDs.min()

            # fixing the wrapping boundaries
            side1_IDs = plate_IDs[0, :]
            side2_IDs = plate_IDs[-1, :]

            for i in range(plate_IDs.shape[1]):
                if side1_IDs[i] != side2_IDs[i]:
                    plate_IDs[plate_IDs == side2_IDs[i]] = side1_IDs[i]

            raw_unified_IDs = np.unique(plate_IDs)
            for i, ID in enumerate(raw_unified_IDs):
                plate_IDs[plate_IDs == ID] = -i
            plate_IDs *= -1
            
        elif unification_axis == 1:
            # fixing the polar condition
            pole1_IDs = np.unique(plate_IDs[0, :])
            for ID in pole1_IDs:
                plate_IDs[plate_IDs == ID] = pole1_IDs.min()
            pole2_IDs = np.unique(plate_IDs[-1, :])
            for ID in pole2_IDs:
                plate_IDs[plate_IDs == ID] = pole2_IDs.min()

            # fixing the wrapping boundaries
            side1_IDs = plate_IDs[:, 0]
            side2_IDs = plate_IDs[:, -1]

            for i in range(plate_IDs.shape[0]):
                if side1_IDs[i] != side2_IDs[i]:
                    plate_IDs[plate_IDs == side2_IDs[i]] = side1_IDs[i]

            raw_unified_IDs = np.unique(plate_IDs)
            for i, ID in enumerate(raw_unified_IDs):
                plate_IDs[plate_IDs == ID] = -i
            plate_IDs *= -1
        else:
            raise ValueError("unification_axis must be either 0 or 1.")
            

    return plate_IDs
        