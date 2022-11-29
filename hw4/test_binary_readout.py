# ========================================
#  QIC  >>  assignment 4
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 23 November 2022
# ========================================

from scipy.io import FortranFile
import numpy as np

#  This script is to test the fortran matrix readout

def read_real8_matrix(fname : str) -> np.matrix:
    """Read binary real8 matrix file"""
    f = FortranFile( fname, 'r')
    shape = f.read_ints(np.int32)
    print( "reading ({},{}) matrix".format(*shape) )
    mx = f.read_reals(float)
    return mx.reshape(shape, order="F")

lol = read_real8_matrix('data/1000-10.0-1.0.dat')

# print head(3) of matrix
print( lol[0:2,0:2] )