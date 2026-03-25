![logo](https://github.com/pjavaheri/platerecipy/raw/master/logo.png)

# `platerecipy`: a package for PLATE RECognition In PYthon
`platerecipy` is a tool for detecting candidate plates on global geophysical 
datasets. It analyzes the surface to identify diffuse and non-conforming regions, 
as well as regions with low confidence in plate assignment.

## Supported platforms

Though `platerecipy` is fundamentally a `Python` package, it also relies on 
backend functionalities that are implemented in `C`. **Linux** is the main operating 
system for `platerecipy`. However, expect for a number of optimizations (e.g., multithreaded `C` functions), 
`platerecipy` is now available via `pip` on Windows systems.

## Installation

Releases are made available on PyPI so that `platerecipy` can be installed using 
`pip` on the shell as follows:
```bash
python -m pip install platerecipy
```

## User guide

`platerecipy`'s documentation is moved to [platerecipy.readthedocs.io](https://platerecipy.readthedocs.io/)!

For a detailed explanation and illustration of `platerecipy`'s core functionalities,
please refer to:
> Javaheri, P., & Lowman, J. P. (2026). A random walker algorithm for plate boundary detection in spherical mantle convection models and global geophysical data sets: Application to Euler vector determination. _Journal of Geophysical Research: Solid Earth_, 131, e2025JB032259. https://doi.org/10.1029/2025JB032259

As a demonstration, assuming a given `input_xs`, `input_ys`, `input_zs`, and  
`input_field` (all `numpy` arrays), using `platerecipy` is as simple as follows:
```python
from platerecipy.model import PlateModel
from platerecipy.grid import SphericalGrid

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

# outputting as a ParaView readable .vtk file
from platerecipy import io
io.save_as_vtk(m)
```

--------------------------------------------------------------------------------

Developed by Pejvak Javaheri, 
[pejvak.javaheri@mail.utoronto.ca](pejvak.javaheri@mail.utoronto.ca).


