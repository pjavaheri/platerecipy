![Version](https://img.shields.io/github/v/tag/pjavaheri/platerecipy?label=Version&color=blueviolet)
[![Version](https://img.shields.io/github/v/tag/pjavaheri/platerecipy?label=PyPI&color=red)](https://pypi.org/project/platerecipy/)
[![https://img.shields.io/badge/License-MIT-yellow.svg](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/pjavaheri/platerecipy/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/Documentation-readthedocs-blue.svg)](https://platerecipy.readthedocs.io/)
![Integration Tests](https://img.shields.io/github/actions/workflow/status/pjavaheri/platerecipy/tests.yml?branch=master&label=Integration%20Tests)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15625715.svg)](https://doi.org/10.5281/zenodo.15625715)

![logo](https://github.com/pjavaheri/platerecipy/raw/master/logo.png)

# `platerecipy`: a package for PLATE RECognition In PYthon
`platerecipy` is a tool for detecting candidate plates on global geophysical 
datasets. It analyzes the surface to identify diffuse and non-conforming regions, 
as well as regions with low confidence in plate assignment.

## Supported platforms

`platerecipy` is now fully available on Linux, Windows, and MacOS. Yet, Linux remains 
the main operating system for testing and official support due to a number of backend functionalities implemented in `C`. 

## Installation

Version releases are made available on PyPI so that `platerecipy` can be installed using 
`pip` on the shell as follows:
```bash
pip install platerecipy
```

## Quick start

`platerecipy`'s documentation is moved to [platerecipy.readthedocs.io](https://platerecipy.readthedocs.io/)!

For a detailed explanation of `platerecipy`'s core methodology,
please refer to:
> Javaheri, P., & Lowman, J. P. (2026). A random walker algorithm for plate boundary detection in spherical mantle convection models and global geophysical data sets: Application to Euler vector determination. _Journal of Geophysical Research: Solid Earth_, 131, e2025JB032259. https://doi.org/10.1029/2025JB032259

**`platerecipy` is capable of operating on both ungridded (e.g., geophysical measurments) and gridded (e.g., geodynamic model output) data sets. Ungridded data are automatically interpolated to a grid of choice, while for a gridded input, plate detection takes places directly on the input grid struture.**

### Ungridded data
As a demonstration, assuming a given `input_xs`, `input_ys`, `input_zs`, and  
`input_field` (all `numpy` arrays), using `platerecipy` is as simple as follows:
```python
from platerecipy.grid import SphericalGrid
from platerecipy.model import PlateModel

input_xs    = # to be specified ...
input_ys    = # to be specified ...
input_zs    = # to be specified ...
input_field = # to be specified ...

# generating a consistent grid for interpolation
grid = SphericalGrid(input_xs, input_ys, input_zs)

# interpolating an input field
field = grid.interpolate_field(input_field)

# initializing a plate model
m = PlateModel(grid)

# stacking the interpolated field
m.stack_field(field, take_log=True)

# finding plates on the stacked field
m.find_plates(
    boundary_quantile     = 0.9,            # threshold for the boundaries 
    separation_tolerance  = 4*3.1416/180.,  # 4 degrees for separation tolerance
    RW_beta               = 200,            # RW beta (for feature sharpness)
    min_marker_size       = 100             # to filter out micro plates
)

m.plate_IDs   # --> the array of plate IDs
m.ID_probs    # --> the array of ID probabilities
```

### Gridded data
If the input data set is already gridded (such as model output), the grid can be directly imported
to bypass interpolation, thus increasing processing speed and accuracy. To do so, the data set must 
be contained in a `.vtu` file and read using the `CustomGrid` object as follows:

```python
from platerecipy.grid import CustomGrid
from platerecipy.model import PlateModel

grid = CustomGrid("dataset.vtu")

# initializing a plate model
m = PlateModel(grid)

# stacking a sample field (must be present in the same vtu file)
m.stack_field(grid.mesh["Strain-rate"], take_log=True)

# finding plates on the stacked field
m.find_plates(
  ...
)
```



## Citing `platerecipy`
If this package has been useful to your research, please cite the original paper:

> Javaheri, P., & Lowman, J. P. (2026). A random walker algorithm for plate boundary detection in spherical mantle convection models and global geophysical data sets: Application to Euler vector determination. _Journal of Geophysical Research: Solid Earth_, 131, e2025JB032259. https://doi.org/10.1029/2025JB032259
```
@article{platerecipy,
author = {Javaheri, P. and Lowman, J. P.},
title = {{A random walker algorithm for plate boundary detection in spherical mantle convection models and global geophysical data sets: Application to Euler vector determination}},
journal = {Journal of Geophysical Research: Solid Earth},
volume = {131},
number = {3},
pages = {e2025JB032259},
keywords = {plate tectonics, plate modeling, Euler vector, diffuse zones, mantle convection, strain-rate},
doi = {https://doi.org/10.1029/2025JB032259},
url = {https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2025JB032259},
eprint = {https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2025JB032259},
note = {e2025JB032259 2025JB032259},
year = {2026}
}
```

as well as the appropriate software version DOI released on `Zenodo`: 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15625715.svg)](https://doi.org/10.5281/zenodo.15625715)

--------------------------------------------------------------------------------

Developed by Pejvak Javaheri, 
[pejvak.javaheri@mail.utoronto.ca](pejvak.javaheri@mail.utoronto.ca).


