# ========================================
#  QIC  >>  assignment 6
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 9 December 2022
# ========================================



# %% imports

import numpy as np
from qiskit import quantum_info as qi
import os
from scipy.io import FortranFile

query_prefix = 'data/test/'
answs_prefix = 'data/answ/'

tolerance = 1e-5


def read_density_matrix(fname) -> dict:
    # example:  read_density_matrix('data/test/lol.dat')
    
    meta = {}

    f = FortranFile( fname, 'r')
    dim = f.read_ints(np.int32)[0]
    nn  = f.read_ints(np.int32)[0]
    meta['dim'], meta['nn'] = dim, nn

    #print( " [file] reading structure (dim = {}, nn = {})".format(dim, nn) )

    #mx = np.empty( (dim**nn,dim**nn), dtype=np.cdouble)
    mx = f.read_record( np.cdouble )
    
    meta['data'] = np.reshape( mx, (dim**nn,dim**nn) )
    meta['label'] = ''
    return meta


# %% execute the test

# look for files in answer directory
data_files = [f for f in os.listdir(answs_prefix) if f.endswith('.dat')]
print( "found {} records in {}".format(len(data_files), answs_prefix) )


for dfile in sorted(data_files):

    label = dfile.split("+", 1)
    if len(label) == 1:   trace, label = 1, label[0]
    else:                 trace, label = int( label[0] ) - 1, label[1]
    print( "loading {}".format(dfile))

    # load the matrix to test
    test_meta = read_density_matrix( query_prefix + label)
    
    # compute the partial trace with qiskit
    rho = qi.partial_trace( test_meta['data'],
                             [ i for i in range(test_meta['nn']) if i != trace] )

    # load the partial trace computed with FORTRAN
    answ_meta = read_density_matrix( answs_prefix + dfile)

    test_result = np.allclose(rho, answ_meta['data'] , rtol=1e-04)
    if test_result:
        print("SUCCESSFUL TEST")
    else:
        print("FAILED TEST")

    print('')

# %%
