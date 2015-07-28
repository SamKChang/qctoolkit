from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

c_module = [Extension(name = 'qctoolkit.demo',
              sources = ['qctoolkit/src/demo.c']),
            Extension(name = "qctoolkit.read_cube", 
              sources = ['qctoolkit/src/readcubemodule.c']),
            Extension(name = 'qctoolkit.primesc',
              sources = ['qctoolkit/src/primesc.pyx'])
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
    'qctoolkit'
  ],
  ext_modules = cythonize(c_module)
)
