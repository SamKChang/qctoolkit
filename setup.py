from setuptools import setup
from distutils.core import Extension

module1 = Extension('demo',
                    sources = ['qctoolkit/src/demo.c'])

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
  ext_modules = [module1]
)
