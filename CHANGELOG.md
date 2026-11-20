# Changelog

Summary of key changes from the backend to the user interface.

## Version 2.x

### Version 2.2.4
* Minor bug fix

### Version 2.2.3
* Updating the order of operations in PlateModel for consistency

### Version 2.2.2
* Fixing C-type pointer declerations

### Version 2.2.1
* Patching function exports (DLL) from grid.h on Windows

### Version 2.2.0
* Embedding version number and segmentation parameters in output `vtu` files
* Added `fast_mode` as default to reuse markers if segmentation parameters are unchanged from a previous call
* Changed default `tangent-plane` interpolation parameters to adjust the stablization value based on the Frobenius norm.
* More examples on the documentaition
* Minor bug fixes 

### Version 2.1.0
* Direct import of irregular grids via `vtu` files
* Automatic surface extraction of spherical shells
* A comprehensive set of integration tests added
* Automated GitHub workflow for running tests
* Minor bug fixes

### Version 2.0.5
* Fixing Github workflow issues

### Version 2.0.4
* Fixing Github workflow issues

### Version 2.0.3
* Rewriting the backend to store nodal connections explicitly for a more robust and generalized implementation of the Random Walker algorithm
* Using `vtp` files and the `pyvista` package for custom (irregular) spherical grids, primarily to simplify `ASPECT` users workflow
* Using `OpenMP` for parallelization in Linux and Windows, and leaving the `pthread.h` header only for MacOS
* Minor bug fixes

### Version 2.0.2
* Minor bug fix for automatic documentation on readthedocs

### Version 2.0.1
* A number of fixes for Windows compatibility 
* PyPI automatic compilation for Linux, MacOS, and Windows

### Version 2.0.0
* Added partial spherical grid
* Implemented a built-in `vtk` output functionality 
* Changed `pyvista` and `vtk` packages from required to optional dependencies