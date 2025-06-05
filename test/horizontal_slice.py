"""
@file horizontal_slice.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief A script for testing and visualizing the sensitivity to RW connection weights.
"""

# ~~~~~~~~~~~~~~ Reading and preparing an input dataset ~~~~~~~~~~~~~~~~~
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

input_field = np.load('../test/field.npy')
R = 1.
thetas, phis = np.meshgrid(
    np.linspace(0, np.pi, input_field.shape[0]),
    np.linspace(-np.pi, np.pi, input_field.shape[1]),
    indexing='ij'
)
input_xs = R*np.sin(thetas)*np.cos(phis)
input_ys = R*np.sin(thetas)*np.sin(phis)
input_zs = R*np.cos(thetas)

# ~~~~~~ Stacking the input field and initializing a plate model ~~~~~~
from platerecipy.model import PlateModel
from platerecipy.grid import SphericalGrid

# generating a consistent grid for interpolation
grid = SphericalGrid(input_xs, input_ys, input_zs)

# interpolating an input field
field = grid.interpolate_field(input_field)

# initializing a plate model
m = PlateModel(grid)

# stacking the interpolated field
m.stack_field(field)


# ~~~~ Plotting output obtained from different values of beta ~~~~
fig, axes = plt.subplots(3, 3, figsize=(12, 10))

i = 150

ax = axes[0][0]
ax.set_title('The stacked field, $\mathcal{F}$')
im = ax.imshow(m.stacked_field)
fig.colorbar(im, shrink=0.5, orientation='horizontal')
ax.axhline(i, label=f'1D slice at i={i}', color='red')
ax.legend()

beta = 0
# finding plates on the stacked field
m.find_plates(
    boundary_quantile     = 0.8,              # threshold for the boundaries 
    separation_tolerance  = 2.5*3.1416/180.,  # 4 degrees for separation tolerance
    RW_beta               = beta,               # RW beta (for feature sharpness)
    #min_marker_size       = 20,               # to filter out micro plates
    preserve_small_markers= True
)

ax = axes[1][0]
ax.set_title(r'Plate IDs with $\beta=' + f'{beta}$')
bounds = np.linspace(m.plate_IDs.min()-0.5, m.plate_IDs.max()+0.5, m.plate_IDs.max()+1)
norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
im = ax.imshow(m.plate_IDs, norm=norm, cmap='nipy_spectral')
cb = fig.colorbar(im, shrink=0.5, orientation='horizontal')
cb.set_ticks([1, m.plate_IDs.max()])

ax = axes[1][1]
ax.set_title(r'ID probabilities with $\beta=' + f'{beta}$')
im = ax.imshow(m.ID_probs, vmin=0, vmax=1, cmap='coolwarm')
cb = fig.colorbar(im, shrink=0.5, orientation='horizontal')
cb.set_ticks([0, 1])

ax = axes[1][2]
ax.set_title(r'ID probabilities of the 1D slice with $\beta=' + f'{beta}$')
ax.plot(m.stacked_field[i, :], label=r'$\mathcal{F}(\phi)$', color='red')
ax.plot(m.ID_probs[i, :], label=r'$\mathcal{P}(\phi)$', color='green')
ax.legend()

axes[0][1].plot(m.plate_IDs[i, :], label=r'$\beta =' + f'{beta}$')
axes[0][2].plot(m.ID_probs[i, :], label=r'$\beta =' + f'{beta}$')

beta = 50
# finding plates on the stacked field
m.find_plates(
    boundary_quantile     = 0.8,              # threshold for the boundaries 
    separation_tolerance  = 2.5*3.1416/180.,  # 4 degrees for separation tolerance
    RW_beta               = beta,               # RW beta (for feature sharpness)
    #min_marker_size       = 20,               # to filter out micro plates
    preserve_small_markers= True
)

ax = axes[2][0]
ax.set_title(r'Plate IDs with $\beta=' + f'{beta}$')
bounds = np.linspace(m.plate_IDs.min()-0.5, m.plate_IDs.max()+0.5, m.plate_IDs.max()+1)
norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
im = ax.imshow(m.plate_IDs, norm=norm, cmap='nipy_spectral')
cb = fig.colorbar(im, shrink=0.5, orientation='horizontal')
cb.set_ticks([1, m.plate_IDs.max()])

ax = axes[2][1]
ax.set_title(r'ID probabilities with $\beta=' + f'{beta}$')
im = ax.imshow(m.ID_probs, vmin=0, vmax=1, cmap='coolwarm')
cb = fig.colorbar(im, shrink=0.5, orientation='horizontal')
cb.set_ticks([0, 1])

ax = axes[2][2]
ax.set_title(r'ID probabilities of the 1D slice with $\beta=' + f'{beta}$')
ax.plot(m.stacked_field[i, :], label=r'$\mathcal{F}(\phi)$', color='red')
ax.plot(m.ID_probs[i, :], label=r'$\mathcal{P}(\phi)$', color='green')
ax.legend()

ax = axes[0][1]
ax.set_title(r'Plate IDs along the 1D slice')
ax.plot(m.plate_IDs[i, :], label=r'$\beta =' + f'{beta}$')
ax.legend()

ax = axes[0][2]
ax.set_title(r'ID probabilities along the 1D slice')
ax.plot(m.ID_probs[i, :], label=r'$\beta =' + f'{beta}$')
ax.legend()

fig.tight_layout()
fig.savefig('horizontal_slice.png', dpi=300)