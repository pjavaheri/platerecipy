from platerecipy import grid, model, io

# interpolating the data onto a uniform spherical grid
field, coords = grid.interpolate_to_spherical(
    xs, ys, zs,                                 # Cartesian coordinates
    edots[:,:,-1], take_log = True              # the field, changing in orders
)                                               # of magnitude

# creating an empty model
m = model.PlateModel()

# stacking the interpolated field
m.stack_field(field, take_log=True)

# detecting plates
m.find_plates(
    interior_quantile       = 0.6, 
    wraparound_azimuthally  = True, 
    min_marker_size         = 100
)

# storing plate DIs in various formats
io.save_as_vtp(m, other_fields=[coords[1], coords[2]], other_fields_names=['lat', 'lon'])
io.save_as_csv(m, other_fields=[coords[1], coords[2]], other_fields_names=['lat', 'lon'])
io.save_as_mat(m, other_fields=[coords[1], coords[2]], other_fields_names=['lat', 'lon'])