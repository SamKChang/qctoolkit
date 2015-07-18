import numpy
import ctypes
from numpy.ctypeslib import ndpointer
import ctypes
import glob, os, re

# general path to c_extension.pyc script
lib_path = re.sub(
             'c_extensions\.pyc', 
             '', 
             os.path.realpath(__file__))

lib_demo = ctypes.cdll.LoadLibrary(lib_path + "../demo.so")
fun_demo = lib_demo.cfun
fun_demo.restype = None
fun_demo.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                     ctypes.c_size_t,
                     ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def democ(indata, outdata):
    assert indata.size == outdata.size
    fun_demo(indata, indata.size, outdata)
