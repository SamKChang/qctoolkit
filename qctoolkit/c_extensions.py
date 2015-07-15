import numpy
import ctypes
from numpy.ctypeslib import ndpointer
import ctypes
import glob, os, re

#print glob.glob('../build/lib*/')
lib_path = re.sub(
             'c\_extensions\.pyc', 
             '', 
             os.path.realpath(__file__))
#lib_path = glob.glob('../build/lib*/')[0]
#demo_path = os.path.join(lib_path, "demo.so")
lib = ctypes.cdll.LoadLibrary(lib_path + "../demo.so")
fun = lib.cfun
fun.restype = None
fun.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_size_t,
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def wrap_fun(indata, outdata):
    assert indata.size == outdata.size
    fun(indata, indata.size, outdata)
