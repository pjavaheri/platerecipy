from setuptools import setup, Extension, find_packages

extra_compile_args=[
    '-O3','-std=c99', '-fPIC', '-DLIBCTOOLS'
]

platerecipy_clib_transform_module = Extension(
    'platerecipy_clib_transform',
    sources = [
        'src/clib/transform.c'
    ],
    include_dirs = ['src/clib'],
    extra_compile_args=extra_compile_args,
)

platerecipy_clib_segmentation_module = Extension(
    'platerecipy_clib_segmentation',
    sources = [
        'src/clib/segmentation.c'
    ],
    include_dirs = ['src/clib'],
    extra_compile_args=extra_compile_args,
)

setup(
    name = 'platerecipy',
    version = '1.1.0',
    description = 'PLATE RECognition In PYthon',
    long_description = 
        'platerecipy is a tool for detecting candidate plates on global geophysical ' 
        'datasets. It analyzes the surface to identify diffuse and non-conforming regions, '
        'as well as regions with low confidence in plate assignment.',
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/pjavaheri/platerecipy',
    author = 'Pejvak Javaheri',
    author_email = 'pejvak.javaheri@mail.utoronto.ca',
    license = 'MIT',
    packages = find_packages(where='src', include=['platerecipy']),
    package_dir = {'' : 'src'},
    install_requires = [
        'numpy', 
        'scipy',
        'matplotlib', 
        'pyvista',
        'pandas'
    ],
    ext_modules = [
        platerecipy_clib_transform_module,
        platerecipy_clib_segmentation_module
    ],
    zip_safe = False
)