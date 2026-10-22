"""
@file test_integration_spherical.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Integration test for a full spherical grid.
"""

import numpy as np
from platerecipy.grid import SphericalGrid, CustomGrid
from platerecipy.model import PlateModel
from platerecipy.segmentation import get_AB_from_map
from platerecipy.io import save_as_vtu
import ctypes

## Part 1 (creating a regular grid)

# creating a pseudo dataset 
num_pts = 1000
theta_range = (0.1, np.pi-0.1)
phi_range   = (-np.pi, np.pi)
rng = np.random.default_rng(seed=17)
data = rng.uniform(
    low  =[theta_range[0], phi_range[0]],
    high =[theta_range[1], phi_range[1]],
    size =(num_pts, 2)
)
data_thetas = data[:, 0]
data_phis   = data[:, 1]
data_field  = 3*np.sin(data_thetas)**2 - 1

r = 10.
data_xs = r*np.sin(data_thetas)*np.cos(data_phis)
data_ys = r*np.sin(data_thetas)*np.sin(data_phis)
data_zs = r*np.cos(data_thetas)

# genrating a spherical grid
grid1 = SphericalGrid(data_xs, data_ys, data_zs)
field1 = grid1.interpolate_field(data_field)

# finding plates
model1 = PlateModel(grid1)
model1.stack_field(field1)
model1.find_plates(boundary_absolute = 0.4, RW_beta=100)
save_as_vtu(model1, filename="platerecipy_test_integration_spherical_01.vtu")


## Part 2 (reading the vtu file as a custom grid)

# reading the newly-created vtu file
grid2 = CustomGrid("platerecipy_test_integration_spherical_01.vtu")
field2 = grid2.mesh['stacked_field']

# finding plates
model2 = PlateModel(grid2)
model2.stack_field(field2)
model2.find_plates(boundary_absolute = 0.4, RW_beta=100)
save_as_vtu(model2, filename="platerecipy_test_integration_spherical_02.vtu")

# Part 3 (solver test)
A1, B1, num_labelled1 = get_AB_from_map(
    data             = model1.stacked_field,
    labels           = model1.labels,
    beta             = 100.,
    largest_label    = model1.labels.max(),
    grid             = grid1
)
A2, B2, num_labelled2 = get_AB_from_map(
    data             = model2.stacked_field,
    labels           = model2.labels,
    beta             = 100.,
    largest_label    = model2.labels.max(),
    grid             = grid2
)


## Analysis
print("\n\n\n\n\n\n\n\n ~~~~~~~~~<< Test Begins >>~~~~~~~~~")

print("\n\nChecking if the stacked field is correctly written/read to/from the vtu file ...")
assert np.all(model1.stacked_field.ravel()[field1.shape[1]-1:field1.shape[1]-1+field2.size] == field2), "Expected to be identical. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the labels are the same ...")
assert np.all(model1.labels.ravel()[field1.shape[1]-1:field1.shape[1]-1+field2.size] == model2.labels), "Expected to be identical. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the plate_IDs are correctly reproduced ...")
assert np.all(model1.plate_IDs.ravel()[field1.shape[1]-1:field1.shape[1]-1+field2.size] == model2.plate_IDs), "Expected to be identical. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the ID_probs are correctly reproduced ...")
assert np.all(model1.ID_probs.ravel()[field1.shape[1]-1:field1.shape[1]-1+field2.size] == model2.ID_probs), "Expected to be identical. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the sparse linear system is correctly reproduced ...")
assert np.all(B1 == B2), "Expected RHS vectors to be identical. Integration test failed (x)"
print(" ... RHS vectors are identical")
assert (A1 != A2).nnz == 0, "Expected the LHS sparse matrices to be identical. Integration test failed (x)"
print(" ... LHS sparse matrices are identical")
print("    ... integration test passed!")

print("~~~~~~~~~<< Done >>~~~~~~~~~")