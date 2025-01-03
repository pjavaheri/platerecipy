from platerecipy.model import PlateModel
from platerecipy.grid import SphericalGrid

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
    boundary_quantile  = 0.9,          # threshold for the boundaries 
    # ...
)

# outputting as a ParaView readable .vtp file
from platerecipy import io
io.save_as_vtp(m)