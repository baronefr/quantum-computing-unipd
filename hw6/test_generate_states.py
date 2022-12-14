# ========================================
#  QIC  >>  assignment 6
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 9 December 2022
# ========================================





# %% imports

import numpy as np
from scipy.io import FortranFile
from qiskit import quantum_info as qi

# %% ======================================================
#    I/O definitions
# =========================================================

def write_density_matrix(dim, nn, rho, fname):
    """This function writes a density matrix in a binary file, supported
    by FORTRAN script."""
    f = FortranFile( fname, 'w')
    f.write_record( dim )
    f.write_record( nn )
    f.write_record( rho.T.astype(np.cdouble) )  # recall that Fotran transposes...



query_prefix = 'data/test/'
answs_prefix = 'data/answ/'

fname = 'leel.dat'


query_file = query_prefix + "{}.dat"


# %% ======================================================
#    generate states
# =========================================================


# %% create a simple |0><0|
dim, nn = 2, 2

rho = np.zeros( (dim**nn, dim**nn), dtype=np.cdouble)
rho[0, 0] = 1

write_density_matrix(dim, nn, rho, query_file.format('00') )




# %% create a Bell state |00> + |11>
dim, nn = 2, 2

rho = np.zeros( (dim**nn, dim**nn), dtype=np.cdouble)
rho[0, 0] = 1/2
rho[3, 3] = 1/2
rho[3, 0] = 1/2
rho[0, 3] = 1/2

write_density_matrix(dim, nn, rho, query_file.format('bell0011') )



# %% random density matrix (1)
dim, nn = 2, 5

rho = qi.random_density_matrix( dim**nn )
write_density_matrix(dim, nn, np.matrix( rho , dtype=np.cdouble),
                     query_file.format('rand5') )



# %% random density matrix (2)
dim, nn = 2, 10                # WARNING: nn=14 is \sim 4GB

rho = qi.random_density_matrix( dim**nn )
write_density_matrix(dim, nn, np.matrix( rho , dtype=np.cdouble),
                     query_file.format('rand10') )

# refs -------------
# https://qiskit.org/documentation/stubs/qiskit.quantum_info.random_density_matrix.html
# https://qiskit.org/documentation/stubs/qiskit.quantum_info.partial_trace.html


# %%
