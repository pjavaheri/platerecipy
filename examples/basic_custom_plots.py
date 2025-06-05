"""
@file basic_custom_plots.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief A template for basic usage of platerecipy in customized plots.
"""

# ~~~~~~~~~~~~~~ Reading and preparing an input dataset ~~~~~~~~~~~~~~~~~
import numpy as np

input_field = np.load('field.npy')
R = 1.
thetas, phis = np.meshgrid(
    np.linspace(0, np.pi, input_field.shape[0]),
    np.linspace(-np.pi, np.pi, input_field.shape[1]),
    indexing='ij'
)
input_xs = R*np.sin(thetas)*np.cos(phis)
input_ys = R*np.sin(thetas)*np.sin(phis)
input_zs = R*np.cos(thetas)

# ~~~~~~~~~~ Using platerecipy ~~~~~~~~~~~

from platerecipy.model import PlateModel
from platerecipy.grid import SphericalGrid

input_xs    = input_xs    # to be specified ...
input_ys    = input_ys    # to be specified ...
input_zs    = input_zs    # to be specified ...
input_field = input_field # to be specified ...

# generating a consistent grid for interpolation
grid = SphericalGrid(input_xs, input_ys, input_zs)

# interpolating an input field
field = grid.interpolate_field(input_field, take_log=False)

# initializing a plate model
m = PlateModel(grid)

# stacking the interpolated field
m.stack_field(field, take_log=False)

# finding plates on the stacked field
m.find_plates(
    boundary_quantile      = 0.8,            # threshold for the boundaries 
    separation_tolerance   = 2.5*np.pi/180., # 2 degrees for separation tolerance
                                             # patching 4-degree gaps
    RW_beta                = 50,             # RW beta (for feature sharpness)
    min_marker_size        = 0,              # to potentially filter out micro plates
    preserve_small_markers = True
)

# to create custom plots
import matplotlib.pyplot as plt

# a given plate
ID = 2

fig, ax = plt.subplots()
ax.set_title(
    f"Plate with plate ID = {ID}\n" 
    + f'mean value across the plate = {np.mean(input_field[m.plate_IDs == ID]):7.4f}'
)
ax.imshow(m.plate_IDs == ID)
fig.savefig('basic_custom_plots.png', dpi=300)