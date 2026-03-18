![logo](../media/logo-trans.png)

# Basic Usage

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


