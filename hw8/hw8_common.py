# ========================================
#  QIC  >>  assignment 8
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 24 December 2022
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




# %% ======================================================
#  I/O |  define data structure reading routines
# =========================================================

# data folder
PREFIX = 'data/'
# image folder
img_prefix = "tex/img/"
fname_plots = img_prefix + "{}.pdf"
fname_plots_png = img_prefix + "{}.png"


def read_rgobj(fname : str, records : int = 50) -> dict:
    """Read binary file for RG simulation @ HW8.
    """

    meta = {}

    f = FortranFile( fname, 'r')
    meta['N'] = f.read_ints(np.int64)[0]
    meta['thres'] = f.read_reals(np.float64)[0]
    meta['max_iter'] = f.read_ints(np.int32)[0]

    lambd = []
    counter = []
    gshist = []
    cputime = []
    
    for i in range(records):
        try:
            lambd.append( f.read_reals(np.float64)[0] )
        except:
            print("[info] EOF at record {}".format(i))
            break
        try:
            counter.append( f.read_ints(np.int32)[0] )
            gshist.append( f.read_reals(float) )
            cputime.append( f.read_reals(float)[0] )
        except:
            print("error at record {i}, unexpected EOF")
            return None

    meta['lambd'] = np.array( lambd )
    meta['counter'] = np.array( counter )
    meta['hist'] = np.matrix( gshist )
    meta['cputime'] = np.array( cputime )

    conv = []
    for idx, row in enumerate( meta['hist'] ):
        if( meta['counter'][idx] == -1 ):
            conv.append( row[-1] )
        else:
            conv.append( row.squeeze()[:, meta['counter'][idx] - 1 ] )

    meta['energy'] = np.array( conv ).squeeze() 
    return meta



def read_dmrgobj(fname : str, records : int = 50) -> dict:
    """Read binary file for DMRG simulation @ HW8.
    """

    meta = {}

    f = FortranFile( fname, 'r')
    meta['max_m'] = f.read_ints(np.int32)[0]
    meta['max_iter'] = f.read_ints(np.int32)[0]

    lambd = []
    cputime = []
    counter = []
    energy = []
    hist = []
    hist_err = []

    for i in range(records):
        try:
            lambd.append( f.read_reals(np.float64)[0] )
        except:
            print("[info] EOF at record {}".format(i))
            break

        try:
            cputime.append( f.read_reals(float)[0] )
            counter.append( f.read_ints(np.int32)[0] )
            energy.append( f.read_reals(float)[0] )
            hist.append( f.read_reals(float) )
            hist_err.append( f.read_reals(float) )
        except:
            print("error at record {i}, unexpected EOF".format(i))
            return None

    meta['lambd'] = np.array( lambd )
    meta['cputime'] = np.array( cputime )
    meta['counter'] = np.array( counter )
    meta['energy'] = np.array( energy )
    meta['hist'] = np.array( hist ).squeeze() 
    meta['hist_err'] = np.array( hist_err ).squeeze() 
    return meta




# %% ======================================================
#  PLOT | plot functions
# =========================================================

def mean_field_ising(x) -> float:
    """Mean field theory expectation"""
    if( np.abs(x) <= 2 ):   return -1 - (x/2)**2
    else:                   return -np.abs(x)


def plot_spectrum( this : list, labels : list = None, lw : float = 2, cl: list = None):
    # arg handler
    if isinstance( this, dict):  this = [ this ]

    for ii, tt in enumerate(this):
        if labels is None:
            plt.plot( tt['lambd'], tt['energy'],
                      linewidth = lw
                    )
        else:
            plt.plot( tt['lambd'], tt['energy'],
                      linewidth = lw,
                      label=labels[ii],
                      color = cl[ii]
                    )

    xx = np.linspace(0.0, 3.0, num=1000)
    plt.plot(xx, np.vectorize(mean_field_ising)(xx),
     'k', linestyle=':', label='$MF$')
    
    if labels is not None:  plt.legend()

    plt.ylabel("$\\varepsilon_0$")
    plt.xlabel("$\lambda$")
    return plt



def plot_iteration_to_convergence(this : list, labels : list = None, lw: float = 2):
    # arg handler
    if isinstance( this, dict):  this = [ this ]

    for ii, tt in enumerate(this):
        if labels is None:
            plt.plot( tt['lambd'], tt['counter'],
                      linewidth = lw
                    )
        else:
            plt.plot( tt['lambd'], tt['counter'],
                      label=labels[ii], linewidth = lw  
                    )

    if labels is not None: plt.legend()
    plt.xlabel("$\lambda$")
    return plt


def plot_hist( this : dict, idx : list, lw: float = 2, cm: str = 'plasma'):
    # arg handler
    if isinstance(idx, int):  idx = [ idx ]

    cmap = plt.get_cmap(cm)
    for cc, ii in enumerate(idx):
        y = np.abs( np.diff(this['hist'][ii,:]) )
        plt.plot( np.array(y).squeeze(), 
            label = '${}$'.format(np.round(this['lambd'][ii], 2)),
            linewidth = lw,
            color = cmap(1*(cc)/len(idx))
        )
    
    plt.xlabel('iteration')
    plt.ylabel('improvement per iteration')
    plt.yscale('log')
    plt.legend(title='$\lambda$')
    return plt
# %%
