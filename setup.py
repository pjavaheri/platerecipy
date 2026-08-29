from setuptools import setup, Extension, find_packages
import os           # for PyPI long description 
import sys          # for platform
import sysconfig    # for extension

shared_obj_ext  = sysconfig.get_config_var('EXT_SUFFIX')
if shared_obj_ext is None:
    shared_obj_ext = ".so"

extra_link_args = {'grid': [], 'transform': [], 'segmentation': []}
if sys.platform == 'win32':
    # Windows compile arguments
    extra_compile_args=['-D_GNU_SOURCE']
    extra_link_args['transform'] += ['-fopenmp']
    
elif sys.platform == 'darwin':
    # MacOS compile arguments
    config_vars     = sysconfig.get_config_vars()
    config_vars['LDSHARED'] = config_vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args = {
        'grid': ['-Wl,-install_name,@rpath/libplaterecipy_grid' + shared_obj_ext],
        'transform': ['-Wl,-install_name,@rpath/libplaterecipy_transform' + shared_obj_ext],
        'segmentation': ['-Wl,-install_name,@rpath/libplaterecipy_segmentation' + shared_obj_ext]
    }
    extra_compile_args = [
        '-std=c99',
        '-Wno-unknown-pragmas', 
        '-D_GNU_SOURCE', 
        '-fPIC', 
        '-O3'
    ]
else:
    # Linux compile arguments
    extra_compile_args = [
        '-std=c99',
        '-Wno-unknown-pragmas', 
        '-D_GNU_SOURCE', 
        '-fPIC', 
        '-O3',
        '-fopenmp'
    ]
    extra_link_args['transform'] += ['-fopenmp']

libplaterecipy_grid_module = Extension(
    'libplaterecipy_grid',
    sources = ['src/clib/grid.c'],
    include_dirs        = ['src/clib'],
    define_macros       = [('LIBPLATERECIPY_GRID', None)],
    extra_link_args     = extra_link_args['grid'],
    extra_compile_args  = extra_compile_args + ['-DLIBPLATERECIPY_GRID'],
)

libplaterecipy_transform_module = Extension(
    'libplaterecipy_transform',
    sources = ['src/clib/transform.c'],
    include_dirs        = ['src/clib'],
    define_macros       = [('LIBPLATERECIPY_TRANSFORM', None)],
    extra_link_args     = extra_link_args['transform'],
    extra_compile_args  = extra_compile_args + ['-DLIBPLATERECIPY_TRANSFORM'],
)

libplaterecipy_segmentation_module = Extension(
    'libplaterecipy_segmentation',
    sources = ['src/clib/segmentation.c'],
    include_dirs        = ['src/clib'],
    define_macros       = [('LIBPLATERECIPY_SEGMENTATION', None)],
    extra_link_args     = extra_link_args['segmentation'],
    extra_compile_args  = extra_compile_args + ['-DLIBPLATERECIPY_SEGMENTATION'],
)

with open(
    os.path.join(
        os.path.dirname(__file__), 
        'README.md'
    ), 
    encoding='utf-8'
    ) as readme_file:
    long_description = readme_file.read()

setup(
    name                            = 'platerecipy',
    version                         = '2.0.5',
    description                     = 'PLATE RECognition In PYthon',
    long_description                = long_description,
    long_description_content_type   = 'text/markdown',
    url                             = 'https://github.com/pjavaheri/platerecipy',
    author                          = 'Pejvak Javaheri',
    author_email                    = 'pejvak.javaheri@mail.utoronto.ca',
    license                         = 'MIT',
    packages = find_packages(where='src', include=['platerecipy']),
    package_dir = {'' : 'src'},
    install_requires = [
        'numpy', 
        'scipy',
        'matplotlib',
        'pyvista'
    ],
    #extras_require={
    #    'vtp': ['vtk', 'pyvista']
    #},
    ext_modules = [
        libplaterecipy_grid_module,
        libplaterecipy_transform_module,
        libplaterecipy_segmentation_module
    ],
    zip_safe = False
)