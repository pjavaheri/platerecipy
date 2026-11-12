![logo](../media/logo-trans.png)

# Overview
`platerecipy` has three main modules:

1. `grid`: for generating/importing grids and interpolation

2. `model`: for defining a plate model, candidate plate detection, and probability assessment

3. `io`: containing input/output functions

Each module is designed to work essentially independet from another. That is to say, while the main purpose of this package is plate detection, the toolkit itself (e.g., gridding and interpolation) can utilized for other unrelated tasks.

The processing pipline can be visualized as follows:

[`grid`] --> [`model`] --> [`io`]

---------------------------------------------------------

## **1. An appropriate `Grid` type must be chosen from `platerecipy.grid`:**

| `Grid` Type | Properties |
|--|--|
| `PlanarGrid` | Flat rectilinear grid, specified by (`x_range`,`y_range`, opt. [`x_res`, `y_res`]) |
| `PartialSphericalGrid` | Partial rectilinear spherical grid, specified by (`theta_range`,`phi_range`, opt. [`theta_res`, `phi_res`]) |
| `SphericalGrid` | Rectilinear spherical grid, specified by (`theta_range`,`phi_range`) |
| `CustomGrid` | Any arbitrary grid, specified by an input `vtu` file |


For example, if the data is stored on a grid and in a `.vtu` file, `CustomGrid` is the most suitable choice as it directly imports the grid (and thus bypassing grid generation and data interpolation).

```python
from platerecipy.grid import CustomGrid

grid = CustomGrid("path/to/the/data.vtu")
```
In this case, various fields stored in the `.vtu` file are accessible through the `mesh` member:
 
```python
field = grid.mesh['Strain-rate']
```
Alternatively, for a set of ungridded (scattered) data points on a flat plane, `PlanarGrid` can be used which generates a flat rectilinear grid contained on the xy-plane:

```python
data_xs = ...
data_ys = ...
data_field = ...
x_range = (..., ...)
y_range = (..., ...)

from platerecipy.grid import PlanarGrid

# generating a grid
grid = PlanarGrid(data_xs, data_ys, x_range=x_range, y_range=y_range)

# interpolating data
field = grid.interpolate(data_field)

...
```

!!! tip "Interpolation"
    For any `grid` instance (except for `CustomGrid`), the function `grid.interpolate_field()` interpolates an ungridded field values across the `grid` points. (`CustomGrid` class allows for importing a grid structure directly without the need to interpolate field values.)

!!! tip "`Grid` attributes"
    For ease of use, coordinates and related arrays are can be accessed via: `grid.xs`, `grid.ys`, `grid.thetas`, `grid.r`, etc.


## **2. Creating a `PlateModel` and identifying candidate plates**

A `model` instance is created around an existing `grid` object:
```python
from platerecipy.model import PlateModel

model = PlateModel(grid)
```

With interpolated fields (see the previous section), or, in case of a `CustomGrid` fields embedded in the input `vtu` file (that are accessible via `grid.mesh['Sample field']`), fields are stacked as follows:
```python
model.stack_field(SampleField)
...
```
Then, candidate plates are found as follows:
```python
plate_IDs = model.find_plates(
    boundary_quantile = ...,
    RW_beta = ...,
    separation_tolerance = ...,
    preserve_small_markers = True/False,
    ...
)
```
**Once `model.find_plates()` is called**, the following attributes from the most recent call can be accessed through the `model` instance:

| Attribute | Description |
|--|--|
| `model.stacked_field` | The stacked field |
| `model.plate_IDs` | The Plate ID array |
| `model.ID_probs` | The corresponding probability array |
| `model.markers` | The markers implied by the plate detection parameters |

!!! tip "Iterating on and adjusting plate detection parameters"
    The segmentation marker field (`model.markers`) provides a useful resource for a more educated adjustment of plate detection parameters.

## **3. Exporting values and plotting plate IDs using `io` utilities**
In addition to the arrays listed above, it is possible to export the results by writing as a `vtu` file as follows:
```python
from platerecipy.io import save_as_vtu

save_as_vtu(model, filename="sample.vtu")
```

| Name | Description |
|--|--|
| `save_as_vtu` | Exporting the mesh and field values as an Unstructured Grid data file |
| `save_as_csv` | Exporting field values as a `csv` file |
| `save_as_mat` | Exporting arrays as an `Octave/MATLAB` binary file |
| `save_mollweide_projection` | Plotting the stacked, plate ID, and probability fields as Mollweide projections |
| `save_six_view_angles` | Plotting the plate ID and probability fields in six 3D view angles |