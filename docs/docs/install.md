![logo](media/logo-trans.png)

# Installing `platerecipy`


## Supported platforms

**Now fully available on Linux, MacOS, and Windows!**

## Easy installation
Releases are made available on PyPI so that `platerecipy` can be installed using 
`pip` on the shell as follows:
```bash
pip install platerecipy
```

## An isolated working environment
Especially on Linux, it is generally recommended to create an isolated virtual 
environment, retaining the OS `Python` intact:

```bash
python3 -m venv virtual_env          # creating a virtual environment
source virtual_env/bin/activate      # loading the virtual environment

# installing platerecipy
python -m pip install platerecipy
```

which then can be accessed by:

```bash
source virtual_env/bin/activate
```