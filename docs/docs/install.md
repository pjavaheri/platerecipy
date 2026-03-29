![logo](media/logo-trans.png)

# Installing `platerecipy`


## Supported platforms

Though `platerecipy` is fundamentally a `Python` package, it also relies on 
backend functionalities that are implemented in `C`. **Linux** is the main operating 
system for `platerecipy`.

!!! tip "`platerecipy` is now available on Windows"
    Expect for a number of optimizations (e.g., multithreaded `C` functions), 
    `platerecipy` is now available via `pip` on Windows systems.

## Easy installation
Releases are made available on PyPI so that `platerecipy` can be installed using 
`pip` on the shell as follows:
```bash
pip install platerecipy
```

to also use functions that require `pyvista` package (namely, `platerecipy.io.save_as_vtp`
or `platerecipy.io.save_six_view_angles`), you can install the relevant optional
dependencies by:
```bash
pip install "platerecipy[vtp]"
```

!!! note "`platerecipy.io.save_as_vtk` vs. `platerecipy.io.save_as_vtp`"
    Both of these functions generate ParaView readable output.  However, `save_as_vtk`
    does not require third party packages of `vtk` and `pyvista`. Accordingly, 
    by default `save_as_vtk` is preferred.

## An isolated working environment
Especially on Linux, it is generally recommended to create an isolated virtual 
environment, retaining the OS `Python` intact:

```bash
python3 -m venv virtual_env          # creating a virtual environment
source virtual_env/bin/activate      # loading the virtual environment

# installing platerecipy and jupyter lab
python -m pip install platerecipy jupyterlab

# or to install the optional dependencies: vtk and pyvista
python -m pip install "platerecipy[vtp]" jupyterlab
```

which then can be accessed by:

```bash
source virtual_env/bin/activate
jupyter lab
```