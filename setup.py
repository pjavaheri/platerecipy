from setuptools import setup, Extension, find_packages

extra_compile_args=[
    '-O3','-std=c99', '-fPIC', '-DLIBCTOOLS'
]

clib_transform_module = Extension(
    'clib_transform',
    sources = [
        'src/clib/transform.c'
    ],
    include_dirs = ['src/clib'],
    extra_compile_args=extra_compile_args,
)

clib_utils_module = Extension(
    'clib_utils',
    sources = [
        'src/clib/utils.c'
    ],
    include_dirs = ['src/clib'],
    extra_compile_args=extra_compile_args,
)

setup(
    name = 'platerecipy',
    version = '0.0.1',
    description = 'PLATE RECognition In PYthon',
    long_description = '',
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/pjavaheri/pyyykit',
    author = 'Pejvak Javaheri',
    author_email = 'pejvak.javaheri@mail.utoronto.ca',
    license = 'MIT',
    packages = find_packages(where='src', include=['pyyykit']),
    package_dir = {'' : 'src'},
    install_requires = [
        'numpy', 
        'scipy',
        'matplotlib', 
        'scikit-image'
    ],
    ext_modules = [
        clib_transform_module, 
        clib_utils_module
    ],
    zip_safe = False
)