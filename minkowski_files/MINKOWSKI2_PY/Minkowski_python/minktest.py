from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import ctypes
import c2raytools as c2t

# the data we want to analyse, a 3D array of floats
uv = c2t.read_cbin('/home/gorgel/Documents/C2ray/dtbox_smth7.96.cbin')

# the size of the data we want to analyse, a 1D array of three integers (same as the -x -y -z options)
usizev = np.array([504,504,504])

# integer: number of bins to use (same as -b option)
numbins = int(10)

#float: lowest value of threshold (same as -l option)
low = float(0)

# float: highest values of threshold (same as -h option)
high = float(10) 

output_array = np.zeros((3,numbins+1,5), dtype=np.float32)

# vsumv - the output, a 3D array of floats. The size is (3,numbins+1,5) 

def mink(input_array, array_size, nr_bins, low, high, output_array):
    lib_file = '/home/gorgel/Dropbox/simon/plugg/masterarbete/minkowski_files/MINKOWSKI2_PY/testpython/Minkowski_python/minkowski_python.so'
    lib = ctypes.cdll.LoadLibrary(lib_file)
    func = lib.minkowski
    return func(ctypes.c_void_p(input_array.ctypes.data), ctypes.c_void_p(array_size.ctypes.data),\
         ctypes.c_int(nr_bins), ctypes.c_float(low) , ctypes.c_float(high) ,ctypes.c_void_p(output_array.ctypes.data))
         
A = mink(uv, usizev, numbins, high, low, output_array)
