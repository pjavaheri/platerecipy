![logo](../media/logo-trans.png)

# Basic Usage

!!! tip "Help function"
    It is recommended to access the most up-to-date information by either calling
    `help(example_function)` or running `?example_function` on a `jupyter` notebook.

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
m.stack_field(field)

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

# Additional examples

<div class="example-grid">
  <a href="../detailed_example/">
    <img src="../../media/field_plots_2.png">
    <h3>Detailed example using spherical harmonics</h3>
  </a>
   <a href="../detailed_example/#low-confidence-and-non-conforming-regions">
    <img src="../../media/lowconf-example.png">
    <h3>Low-confidence regions</h3>
  </a>
   <a href="../detailed_example/#mapping-to-original-ungridded-input-data-points">
    <img src="../../media/original_mapping.png">
    <h3>Mapping to original data points</h3>
  </a>
</div>

