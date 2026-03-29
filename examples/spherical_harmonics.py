import numpy as np

def generate_sample_input(lrange=(1, 4), rng_seed=111, npoints=100_000):
    from scipy.special import sph_harm_y
    from scipy.ndimage import laplace

    thetas, phis = np.meshgrid(
        np.linspace(0, np.pi, 300),
        np.linspace(-np.pi, np.pi, 600),
        indexing='ij'
    )
    rng = np.random.default_rng(rng_seed)
    F = 0
    for l in range(lrange[0], lrange[1]):
        for m in range(-l, l+1):
            F += (rng.uniform()-0.25) * sph_harm_y(l, m, thetas, phis).real
    F = laplace(np.abs(F))
    rand_data = rng.choice(thetas.size, size=npoints, replace=False)
    data_phis = phis.ravel()[rand_data]
    data_thetas = thetas.ravel()[rand_data]
    data_F = F.ravel()[rand_data]
    data_xs = np.sin(data_thetas)*np.cos(data_phis)
    data_ys = np.sin(data_thetas)*np.sin(data_phis)
    data_zs = np.cos(data_thetas)
    return data_xs, data_ys, data_zs, data_F


data_xs, data_ys, data_zs, data_F = generate_sample_input()


from platerecipy.grid import SphericalGrid

# initializing a grid base on input coordinates
grid = SphericalGrid(data_xs, data_ys, data_zs)

# interpolating the input data
field1 = grid.interpolate_field(data_F)

from platerecipy.model import PlateModel

# initializing a plate model
model = PlateModel(grid)

# stacking our interpolated field
model.stack_field(field1)

# additional fields can be stacked here

model.find_plates(
    boundary_quantile    = 0.95, 
    separation_tolerance = 1*np.pi/180.,
    RW_beta              = 100.0,
)


import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(10, 6))

ax = axes[0][0]
ax.set_title("Stacked field")
ax.pcolormesh(grid.phis, np.pi/2 - grid.thetas, model.stacked_field)

ax = axes[1][0]
ax.set_title("Markers")
ax.pcolormesh(grid.phis, np.pi/2 - grid.thetas, model.markers)

ax = axes[0][1]
ax.set_title("Plate IDs")
ax.pcolormesh(grid.phis, np.pi/2 - grid.thetas, model.plate_IDs, cmap="jet")

ax = axes[1][1]
ax.set_title("Probability field")
ax.pcolormesh(grid.phis, np.pi/2 - grid.thetas, model.ID_probs, cmap="coolwarm")

fig.tight_layout()
fig.savefig("field_plots.png", dpi=300)

from platerecipy.io import save_mollweide_projection, save_as_vtk

# generating ParaView readable legacy .vtk
save_as_vtk(model)

# generating .png Mollweide projection 
save_mollweide_projection(model)

# mapping to original input data points
org_IDs = grid.map_to_original_input(model.plate_IDs.ravel())

plt.scatter(grid.original_phis, np.pi/2 - grid.original_thetas, c=org_IDs, marker='.', s=0.1)

plt.savefig("original_mapping.png", dpi=300)