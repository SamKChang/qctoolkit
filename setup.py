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

#            Extension(name = "qctoolkit.read_cube", 
#              extra_compile_args=['-O3'],
#              sources = ['qctoolkit/src/readcubemodule.c']),
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
            Extension(name = "qctoolkit.analysis.read_cube", 
              extra_compile_args=['-O3'],
              sources = ['qctoolkit/analysis/c_extension/'+\
                         'readcubemodule.c']),
            Extension(name = "qctoolkit.ML.coulomb_matrix", 
              extra_compile_args=['-O3'],
              sources = ['qctoolkit/src/coulombmatrixmodule.c',
                         'qctoolkit/src/utilities.c']),
            Extension(name = "qctoolkit.MD.dlist_1", 
              sources = ['qctoolkit/MD/c_extension/'+\
                         'dlist1.c'],
              extra_compile_args=['-fopenmp', '-fpic', '-lm',
                                  '-Wno-write-strings'],
              extra_link_args=['-lgomp', '-shared'],
              include_dirs = [np.get_include()]),
            Extension(name = "qctoolkit.MD.dlist_2", 
              sources = ['qctoolkit/MD/c_extension/'+\
                         'dlist2.c'],
              extra_compile_args=['-fopenmp', '-fpic', '-lm',
                                  '-Wno-write-strings'],
              extra_link_args=['-lgomp', '-shared'],
              include_dirs = [np.get_include()]),
            Extension(name = "qctoolkit.MD.vacf", 
              sources = ['qctoolkit/MD/c_extension/'+\
                         'vacf.c'],
              extra_compile_args=['-fopenmp', '-fpic', '-lm',
                                  '-Wno-write-strings'],
              extra_link_args=['-lgomp', '-shared'],
              include_dirs = [np.get_include()]),
            Extension(name = "qctoolkit.QM.gcint", 
              sources = ['qctoolkit/QM/c_extension/'+\
                         'gcint.c'],
              extra_compile_args=['-fopenmp', '-fpic', '-lm',
                                  '-Wno-write-strings'],
              extra_link_args=['-lgomp', '-shared',
                               '-lgsl', '-lgslcblas'],
              include_dirs = [np.get_include()]),
           ]

setup(name='qctoolkit',
  version='0.1.1',
  description='quantum chemistry tool kit',
  url='https://github.com/SamKChang/qctoolkit.git',
  author='K. Y. S. Chang',
  author_email='ky.sam.chang@gmail.com',
  packages=[
    'qctoolkit.analysis',
    'qctoolkit.utilities',
    'qctoolkit.projects',
    'qctoolkit.projects.Basel',
    'qctoolkit.projects.Basel.p01_AlGaAs',
    'qctoolkit.QM',
    'qctoolkit.QM.qmcode',
    'qctoolkit.QM.tools',
    'qctoolkit.data',
    'qctoolkit.data.elements',
    'qctoolkit.data.PP',
    'qctoolkit.data.PP.cpmd',
    'qctoolkit.ccs',
    'qctoolkit.optimization',
    'qctoolkit.ML',
    'qctoolkit.MD',
    'qctoolkit.MD.mdcode',
    'qctoolkit.MD.trajectory',
    'qctoolkit.alchemy',
    'qctoolkit.properties',
    'qctoolkit'
  ],
  package_data={'': ['elements/elements.yml', 
                     'data/PP/cpmd/*.psp']},
  include_package_data=True,
  ext_modules = cythonize(c_module)
)
