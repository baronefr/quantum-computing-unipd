# ========================================
#  QIC  >>  assignment 7
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 15 December 2022
# ========================================



# %% import libraries ---------------------------------------

import numpy as np
import matplotlib.pyplot as plt

from scipy.io import FortranFile



# %% set the plotting env ---------------------------------

useLaTeX = True  # if True, use LaTeX backend
if useLaTeX:
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)

plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (9,6)

# image folder
img_prefix = "tex/img/"
fname_plots = img_prefix + "{}.pdf"
fname_plots_png = img_prefix + "{}.png"




# %% ======================================================
#  I/O |  define data structure
# =========================================================

# data folder
PREFIX = 'data/'

def read_fortran_matrix(fname : str) -> dict:
    """Read binary file with Hamiltonian matrix for HW7.
    """

    meta = {}

    f = FortranFile( fname, 'r')
    N = f.read_ints(np.int64)[0]
    steps = f.read_ints(np.int32)[0]
    neig = f.read_ints(np.int32)[0]
    
    mx = []
    for i in range(steps):
        try:
            mx.append( f.read_reals(float) )
        except:
            print('err at record', i)
            return None

    meta['N'] = N
    meta['steps'] = steps
    meta['neig'] = neig
    meta['data'] = np.matrix(mx)
    return meta




# %% ======================================================
#  I/O | load the dataset
# =========================================================


available_data = [4,8,12,14,15]
dataset = {}

for n in available_data:
    dataset[str(n)] = read_fortran_matrix('data/{}-30-10.hw7'.format(n))


# %% ======================================================
#  PLOT | plot the spectrum energy of ising model
# =========================================================

def mean_field_ising(x) -> float:
    """theoretical mean field prediction of quantum ising ground state"""
    if( np.abs(x) <= 2 ):   return -1 - (x/2)**2
    else:                   return -np.abs(x)


def plot_spectrum( this : dict, n : int, legend : bool = True, figsize = (10,5)):
    if n > this['neig']:
        print('n is higher than computed eigenvalue')
        return None

    plt.figure( figsize=figsize )

    cmap = plt.get_cmap('plasma')
    for i in range(n):
        plt.plot( np.linspace(0.0, 3.0, num=this['steps']), 
                  this['data'][:,i]/(this['N']), label = "${}$".format(i),
                  linestyle = '-' if(i%2 == 1) else '--',
                  color = cmap(1*(i)/n), marker='x', markersize=4 )

    xx = np.linspace(0.0, 3.0, num=1000)
    plt.plot(xx, np.vectorize(mean_field_ising)(xx),
             'k', linestyle=':', label='$MF$')

    if legend: plt.legend(title='$\\varepsilon$', ncol=1, loc = 'lower left')

    plt.title( '$N = {}$'.format( this['N'] ) )
    plt.xlabel("$\lambda$")
    plt.tight_layout()
    return plt



# %% make the plots

neig = 5
plot_spectrum( dataset['4'],  neig,
               legend = True, figsize = (6,8)).savefig( fname_plots.format('spec04'), transparent=True)
plot_spectrum( dataset['8'],  neig, legend = False).savefig( fname_plots.format('spec08'), transparent=True)
plot_spectrum( dataset['15'], neig, legend = False).savefig( fname_plots.format('spec15'), transparent=True) # replace with 15






# %% ======================================================
#  PLOT | plot the energy gap of ising model for various N
# =========================================================

def plot_gap( collection : dict, n : int, low : int = 0, figsize = (8,5) ):
    plt.figure( figsize=figsize )

    cmap = plt.get_cmap('plasma')
    cmapt = len(collection) + 1

    for i, sim in enumerate(collection):
        plt.plot( np.linspace(0.0, 3.0, num=sim['steps']), 
                  sim['data'][:,n] - sim['data'][:,low],
                  label = "${}$".format(sim['N']),
                  color = cmap(1*(i+1)/cmapt),
                  linewidth = 2 )

    if((n == 1) and (low == 0)):
        xx = np.linspace(1.0, 3.0, num=1000)
        plt.plot( xx, 2*(xx - 1),
                  'k', linewidth = 2, 
                  linestyle='--', label='$MF$')
    elif((n == 2) and (low == 0)):
        xx = np.linspace(0.0, 3.0, num=1000)
        plt.plot( xx, np.abs( 2*(xx - 1) ), 
                 'k', linewidth = 2, 
                 linestyle='--', label='$MF$')

    plt.ylabel("$\\varepsilon_{} - \\varepsilon_{}$".format(n, low) )
    plt.xlabel("$\lambda$")
    plt.legend(title='$N$', ncol=2)
    plt.tight_layout()
    return plt



# %% make the plots

# first eigenvalue
plot_gap( dataset.values(), 1, figsize=(6,8)).savefig( fname_plots.format('gap1'), transparent=True)  
# second eigenvalue
plot_gap( dataset.values(), 2, low=1).savefig( fname_plots.format('gap2'), transparent=True)  

# %%
