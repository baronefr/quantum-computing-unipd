# ========================================
#  QIC  >>  assignment 6
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 9 December 2022
# ========================================





# %% imports
import numpy as np
import matplotlib.pyplot as plt

fname_sep = 'data/separable.csv'
fname_gen = 'data/general.csv'
# file format:  dim, nn, time

img_prefix = "tex/img/"
fname_plots = img_prefix + "{}.pdf"


# %% set the plotting env ---------------------------------

useLaTeX = True  # if True, use LaTeX backend
if useLaTeX:
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)

plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (9,6)



# %% ======================================================
#  I/O
# =========================================================


def read_benchmarks(fname : str, data_scale : str) -> dict:
    meta = {}
    dat = np.loadtxt(fname, delimiter=',')

    meta['raw'] = dat
    
    meta['dim'] = np.unique( dat[:,0].astype(int) )
    meta['nn'] = np.unique( dat[:,1].astype(int) )

    meta['time'] = dat[:,2].reshape( len(meta['dim']), len(meta['nn']) )

    if data_scale == '*':
        meta['GB'] = ( (dat[:,0] * dat[:,1]) * 16 / ((1024)**3) ).reshape( len(meta['dim']), len(meta['nn']) )
    elif data_scale == '**':
        meta['GB'] = ( (dat[:,0] ** dat[:,1]) * 16 / ((1024)**3) ).reshape( len(meta['dim']), len(meta['nn']) )
    else:
        print('unhandled data scaling factor')
        meta['GB'] = None
    return meta


# read data from benchmark files
sep = read_benchmarks(fname_sep, data_scale = '*')
gen = read_benchmarks(fname_gen, data_scale = '**')




# %% define the plot functions

def plot_time_mat(metad):
    plt.matshow( metad['time'].T )
    ax = plt.gca()
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('time [$s$]', rotation=270)
    cbar.ax.get_yaxis().labelpad = 24

    plt.ylabel('$N$', rotation=0, labelpad = 15)
    ax.set_yticks( np.arange(len(metad['nn'])) )
    ax.set_yticklabels( ['${}$'.format(i) for i in metad['nn']] )
    plt.xlabel('$D$')
    ax.set_xticks( np.arange(len(metad['dim'])) )
    ax.set_xticklabels( ['${}$'.format(i) for i in metad['dim']] )

    ax.xaxis.set_ticks_position('bottom')
    return plt


def plot_asize(metad : dict, plot_legend : bool = True, linestyle = '-'):
    cmap = plt.get_cmap('inferno')

    for i in range(4,14,2):
        plt.plot(metad['nn'], metad['GB'][i,:], color = cmap(0.8*i/14),
                 linestyle=linestyle,
                 label='${}$'.format(metad['dim'][i]) )
    if plot_legend:
        plt.legend(title = '$D$')
    plt.ylabel('$\\textnormal{allocated memory} \; [GB]$')
    plt.xlabel('$N$')
    plt.yscale('log')
    return plt





# %% ======================================================
#  PLOT |  matrix with CPU time
# =========================================================

plot_time_mat(sep).savefig( fname_plots.format("cpu_sep"), transparent=True)
plot_time_mat(gen).savefig( fname_plots.format("cpu_gen"), transparent=True)



# %% ======================================================
#  PLOT |  allocated GB vs N for increasing D
# =========================================================

plot_asize(gen)
plot_asize(sep, plot_legend = False, linestyle = '--')
plt.savefig( fname_plots.format("malloc_scaling"), transparent=True)


# %%
