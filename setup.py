from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize
import numpy as np

# to enable openmp, use:
#  extra_compile_args=['-fopenmp', '-fpic']
#  and extra_link_args=['-lgomp', '-shared']
# Cython extension:
#  Extension(name = 'qctoolkit.path.to.extension.name',
#    sources = ['qctoolkit/path/to/extension/name.pyx'])

c_module = [Extension(name = "qctoolkit.ML.kernel_matrix", 
              sources = ['qctoolkit/ML/c_extension/'+\
                         'kernelmatrixmodule.c',
                         'qctoolkit/ML/c_extension/kernels.c'],
              extra_compile_args=['-fopenmp', '-fpic',
                                  '-Wno-write-strings'],
              extra_link_args=['-lgomp', '-shared'],
              include_dirs = [np.get_include()]),
            Extension(name = "qctoolkit.ML.kernel_vectors", 
              sources = ['qctoolkit/ML/c_extension/'+\
                         'kernelvectorsmodule.c',
                         'qctoolkit/ML/c_extension/kernels.c'],
              extra_compile_args=['-fopenmp', '-fpic',
                                  '-Wno-write-strings'],
              extra_link_args=['-lgomp', '-shared'],
              include_dirs = [np.get_include()]),
            Extension(name = "qctoolkit.read_cube", 
              extra_compile_args=['-O3'],
              sources = ['qctoolkit/src/readcubemodule.c']),
            Extension(name = "qctoolkit.coulomb_matrix", 
              extra_compile_args=['-O3'],
              sources = ['qctoolkit/src/coulombmatrixmodule.c',
                         'qctoolkit/src/utilities.c'])
           ]

setup(name='qctoolkit',
  version='0.1',
  description='quantum chemistry tool kit',
  url='https://github.com/SamKChang/qctoolkit.git',
  author='K. Y. S. Chang',
  author_email='ky.sam.chang@gmail.com',
  packages=[
    'qctoolkit.projects',
    'qctoolkit.projects.Basel',
    'qctoolkit.projects.Basel.p01_AlGaAs',
    'qctoolkit.io_format',
    'qctoolkit.elements',
    'qctoolkit.ccs',
    'qctoolkit.optimization',
    'qctoolkit.ML',
    'qctoolkit.MD',
    'qctoolkit.alchemy',
    'qctoolkit.properties',
    'qctoolkit'
  ],
  package_data={'': ['elements/elements.yml', 
                     'data/PP/cpmd/*.psp']},
  include_package_data=True,
  ext_modules = cythonize(c_module)
)
