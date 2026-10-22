"""
@file test_integration_partial_spherical.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Integration test for a partial spherical grid.
"""

import numpy as np
from platerecipy.grid import PartialSphericalGrid, CustomGrid
from platerecipy.model import PlateModel
from platerecipy.segmentation import get_AB_from_map
from platerecipy.io import save_as_vtu
import ctypes

## Part 1 (creating a regular grid)

# creating a pseudo dataset 
num_pts = 1000
theta_range = (5*np.pi/180, 15*np.pi/180)
phi_range   = (-90*np.pi/180, 90*np.pi/180)
rng = np.random.default_rng(seed=17)
data = rng.uniform(
    low  =[theta_range[0], phi_range[0]],
    high =[theta_range[1], phi_range[1]],
    size =(num_pts, 2)
)
data_thetas = data[:, 0]
data_phis   = data[:, 1]
data_field  = np.cos(data_phis*10)*data_thetas

r = 10.
data_xs = r*np.sin(data_thetas)*np.cos(data_phis)
data_ys = r*np.sin(data_thetas)*np.sin(data_phis)
data_zs = r*np.cos(data_thetas)

# genrating a spherical grid
grid1 = PartialSphericalGrid(data_xs, data_ys, data_zs, theta_range=theta_range, phi_range=phi_range)
field1 = grid1.interpolate_field(data_field)

# finding plates
model1 = PlateModel(grid1)
model1.stack_field(field1)
model1.find_plates(boundary_absolute = 0.3, RW_beta=100)
save_as_vtu(model1, filename="platerecipy_test_integration_partial_spherical.vtu")


## Part 2 (reading the vtu file as a custom grid)

# reading the newly-created vtu file
grid2 = CustomGrid("platerecipy_test_integration_partial_spherical.vtu")
field2 = grid2.mesh['stacked_field']

# finding plates
model2 = PlateModel(grid2)
model2.stack_field(field2)
model2.find_plates(boundary_absolute = 0.3, RW_beta=100)


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
print(model1.stacked_field-field2.reshape(field1.shape))
assert np.all(model1.stacked_field.ravel() == field2), "Expected all zeros above ^^^. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the labels are the same ...")
print(model1.labels-model2.labels.reshape(field1.shape))
assert np.all(model1.labels == model2.labels.reshape(field1.shape)), "Expected all zeros above ^^^. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the plate_IDs are correctly reproduced ...")
print(model1.plate_IDs-model2.plate_IDs.reshape(field1.shape))
assert np.all(model1.plate_IDs == model2.plate_IDs.reshape(field1.shape)), "Expected all zeros above ^^^. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the ID_probs are correctly reproduced ...")
print(model1.ID_probs-model2.ID_probs.reshape(field1.shape))
assert np.all(model1.ID_probs == model2.ID_probs.reshape(field1.shape)), "Expected all zeros above ^^^. Integration test failed (x)"
print("    ... integration test passed!")

print("\n\nChecking if the sparse linear system is correctly reproduced ...")
assert np.all(B1 == B2), "Expected RHS vectors to be identical. Integration test failed (x)"
print(" ... RHS vectors are identical")
assert (A1 != A2).nnz == 0, "Expected the LHS sparse matrices to be identical. Integration test failed (x)"
print(" ... LHS sparse matrices are identical")
print("    ... integration test passed!")

print("~~~~~~~~~<< Done >>~~~~~~~~~")