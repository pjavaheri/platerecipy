![logo](logo.png)

# `platerecipy`: a package for PLATE RECognition In PYthon

## Supported platforms

Though `platerecipy` is fundamentally a Python package, it also relies on 
backend functionalities that are implemented in C. Typically, Windows compared 
to Linux requires a slightly different compilation and linking which at this 
stage is not supported. That is why the only supported platform is Linux.

## Installation

Until the packages is released on PyPI, `platerecipy` can be installed using 
`pip` on the shell as follows:
```bash
foo@bar:~$ python -m pip install dist/platerecipy-?.?.?.tar.gz
```

```python
import numpy as np
a = np.array([])
```
where `platerecipy-?.?.?.tar.gz` is the installation tar file for the desired 
version. Note that all installation tar files can be found in `dist` 
subdirectory.

## User guide

The user manual can be found under `doc` as a collection of LaTeX files. It 
includes general information about the package, function/class documentations, 
as well as a number of examples to be used as starting recipes! 

--------------------------------------------------------------------------------

Developed by Pejvak Javaheri, 
[pejvak.javaheri@mail.utoronto.ca](pejvak.javaheri@mail.utoronto.ca).


