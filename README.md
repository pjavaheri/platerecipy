![logo](logo.png)

# `platerecipy`: a package for PLATE RECognition In PYthon
`platerecipy` is a tool for detecting candidate plates on global geophysical 
datasets. It analyzes the surface to identify diffuse and non-conforming regions, 
as well as regions with low confidence in plate assignment.

## Supported platforms

Though `platerecipy` is fundamentally a `Python` package, it also relies on 
backend functionalities that are implemented in `C`. Typically, Windows compared 
to Linux requires a slightly different compilation and linking which at this 
stage is not supported. That is why, currently, Linux is the only supported platform.

## Installation

Until the packages is released on PyPI, `platerecipy` can be installed using 
`pip` on the shell as follows:
```bash
user@shell:~$ python -m pip install dist/platerecipy-?.?.?.tar.gz
```
where `platerecipy-?.?.?.tar.gz` is the installation tar file for the desired 
version. Note that all installation tar files can be found in `dist` 
subdirectory.

## User guide

The user manual can be found under `doc` as a LaTeX source file. It 
includes general information about the package, function/class documentations, 
as well as a number of examples to be used as recipes!

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

# outputting as a ParaView readable .vtp file
from platerecipy import io
io.save_as_vtp(m)
```

--------------------------------------------------------------------------------

Developed by Pejvak Javaheri, 
[pejvak.javaheri@mail.utoronto.ca](pejvak.javaheri@mail.utoronto.ca).


