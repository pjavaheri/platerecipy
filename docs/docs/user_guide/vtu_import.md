![logo](../media/logo-trans.png)

# `vtu` Imports and Irregular Grids

`platerecipy` now recognizes irregular grids through `.vtu` files. This is a significant change
since the previous version and allows for an easier incorporation of `platerecipy` 
into workflows that primarily utilize `.vtu` files (e.g., [ASPECT](https://aspect.geodynamics.org/) output).

This functionality is enabled via the `CustomGrid` object that utilizes the `pyvista` package to perform I/O. 

!!! tip "Skipping interpolation"
    With a mesh embedded in the input `vtu` file, `platerecipy` parses faces to generate RW connections. 
    Since the grid is directly read and imported, no interpolation is performed which enhances speed and accuracy.

Considering that the file `sample.vtu` contains
**only the surface of the model** and it also stores the field `Field Name`, 
the following basic template can be used to read ana analyze the data:

```python
from platerecipy.model import PlateModel
from platerecipy.grid import CustomGrid

# reading the input mesh
vtu_file_address = 'sample.vtu'
grid = CustomGrid(vtu_file_address)

# initializing a plate model
m = PlateModel(grid)

# stacking the interpolated field
m.stack_field(grid.mesh['Field Name'])

# finding plates on the stacked field
m.find_plates(
    boundary_quantile     = 0.9,            # threshold for the boundaries 
    separation_tolerance  = 4*3.1416/180.,  # 4 degrees for separation tolerance
    RW_beta               = 200,            # RW beta (for feature sharpness)
    min_marker_size       = 100             # to filter out micro plates
)
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

