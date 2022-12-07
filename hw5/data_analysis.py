# ========================================
#  QIC  >>  assignment 5
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 30 November 2022
# ========================================



# %% import libraries ---------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib as mpl
from scipy.optimize import curve_fit

from scipy.io import FortranFile
import os



# %% set the plotting env ---------------------------------

useLaTeX = True  # if True, use LaTeX backend
if useLaTeX:
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)

plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (9,6)




# %% ======================================================
#  I/O |  define data structure
# =========================================================

# data folder
PREFIX = 'data/'
# image folder
img_prefix = "tex/img/"
fname_plots = img_prefix + "{}.pdf"
fname_plots_png = img_prefix + "{}.png"
do_writeplot = True

def read_fortran_hw5(fname : str):
    """Read binary file for HW5 file structure.

    # FILE STRUCTURE ---------------------
    #
    #  (integer) Nx
    #  (integer) Nt
    #  (real)    xmax
    #  (real)    tmax
    #  (real)    omega
    #  (real)    tau
    #  {  (real) time    (double complex)  psi  }  x Nt
    # ---------------------------------------
    """
    meta = {}

    f = FortranFile( fname, 'r')
    Nx = f.read_ints(np.int32)[0]
    Nt = f.read_ints(np.int32)[0] + 1
    meta['Nx'], meta['Nt'] = Nx, Nt

    meta['xmax'] = f.read_reals(float)[0]
    meta['tmax'] = f.read_reals(float)[0]
    meta['omega'] = f.read_reals(float)[0]
    meta['tau'] = f.read_reals(float)[0]

    print( " [file] reading structure ({},{})".format(Nx,Nt) )

    mx = np.empty( (Nt,Nx), dtype=np.cdouble)
    mt = np.empty(Nt, dtype=float)
    for i in range(Nt):
        mt[i] = f.read_reals(float)[0]
        mx[i,:] = f.read_record( np.cdouble )

    meta['t'] = mt
    return mx, meta


def print_metadata(meta : dict):
    """Prints the metadata associated to a binary HW5 file."""
    print( "[time] :  tmax={} Nt={}".format(meta['tmax'], meta['Nt']) )
    print( "[ x ]  :  xmax={} Nx={}".format(meta['xmax'], meta['Nx']) )
    print( "[hami] :  omega={} tau={}".format(meta['omega'], meta['tau']) )
    




# %% ======================================================
#  DATA | load a specific file for sampling
# =========================================================

data, meta = read_fortran_hw5('data/psi_10001-10000-10.0.dat')
print_metadata(meta)


# %% ======================================================
#  PLOT |  cmap of density as function of time
# =========================================================

def plot_density_timeline(data, meta, cmap='hot',
                          ylim=None, figsize = (14,6), aspect = 0.5):
    # viridis colormap is also good
    plt.figure( figsize = figsize )
    psisq = np.real( data * data.conj() ).T
    plt.imshow( psisq, cmap=cmap,
                interpolation='none', aspect = aspect,
                extent=[0, meta['tmax'], meta['xmax'], -meta['xmax']] )

    if ylim is not None: plt.ylim(*ylim)

    plt.ylabel('$x$', rotation=0)
    plt.xlabel('$t$')
    plt.tight_layout()
    plt.axvline(x = meta['tau'], ymin=0.9, ymax=1, color='r', linestyle='--')
    if ylim is None:
        plt.text(meta['tau'] + 0.2, -0.8*meta['xmax'],
          "$\\tau = {}$".format(round(meta['tau'])), color='w')
    else:
        plt.text(meta['tau'] + 0.2, 0.8*ylim[1],
          "$\\tau = {}$".format(round(meta['tau'])), color='w')
    return plt


def add_fit(data, meta):
    def func(x, m):
        return m*x

    thex = np.linspace(0, meta['tmax'], num=meta['Nt'])
    theline = -meta['xmax'] + np.argmax(
        np.real( data * data.conj() ).T,axis=0
        )*2*meta['xmax']/meta['Nx']

    plt.plot( thex, theline, 'k',
              linewidth=1, linestyle = '--')
    popt, pcov = curve_fit(func, thex, theline)
    plt.plot(thex, func(thex, *popt), 'deeppink',
            label='fit: m=%5.3f' % tuple(popt), linewidth=2 )
    print('fit : m = %5.7f' % tuple(popt), ' +- {}'.format(np.sqrt(np.diag(pcov))))


pplt = plot_density_timeline(data, meta, ylim=(-3,6), aspect=0.8)

if do_writeplot: 
    pplt.savefig( fname_plots_png.format("optimal_density"), transparent=True)




# %% ======================================================
#  PLOT |  residual on linearized drift
# =========================================================

def plot_linearized_residual(data, meta, figsize = (8, 3), matching = True, ylim=None):
    def func(x, m):
        return m*x
    plt.figure(figsize = figsize)
    thex = np.linspace(0, meta['tmax'], num=meta['Nt'])
    theline = -meta['xmax'] + np.argmax(
        np.real( data * data.conj() ).T,axis=0
        )*2*meta['xmax']/meta['Nx']
    popt, pcov = curve_fit(func, thex, theline)

    harm = theline-func(thex, *popt)
    plt.plot(thex, harm, 'k')
    if matching:
        component = np.sin(meta['omega']*thex)
        plt.plot(thex, component*np.max(harm),
                 'r--', linewidth = 1, label="$\sin(\omega t), \omega = {}$".format(int(meta['omega'])) )
        plt.legend(loc='lower right', prop={'size': 18})

    plt.title('$\Psi(t) - \mathcal{L}(\Psi(t))$')
    plt.xlabel('$t$')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    if ylim is not None:
        plt.ylim(ylim)
    plt.tight_layout()
    return plt

pplt = plot_linearized_residual(data, meta)

if do_writeplot: 
    pplt.savefig( fname_plots.format("optimal_density_res"), transparent=True)



# %% ======================================================
#  PLOT |  Real part of wavefunction
# =========================================================

def plot_wave_real(data, meta, tmaxidx, n,
                   basp = 40,
                   figsize = (9, 7), cmapname='plasma'):
    cmap = plt.get_cmap(cmapname)
    fig, ax1 = plt.subplots(1, 1, figsize=figsize)

    timevals = []
    for i in range(0, tmaxidx, int(tmaxidx/n) ):
        ax1.plot(np.linspace(-meta['xmax'], meta['xmax'], num=meta['Nx']),
                 np.real(data[i,:]),
                 color = cmap(i/(tmaxidx-1))
                )
        timevals.append(meta['t'][i])

    norm = mpl.colors.Normalize(vmin=timevals[0], vmax=timevals[-1])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    plt.tight_layout()
    plt.xlabel("$x$", labelpad=2)
    plt.ylabel("$Re(\psi)$")
    plt.colorbar(sm, ticks = timevals, ax = ax1,
                 orientation = 'horizontal',
                 aspect = basp,
                 label = '$t$')
    return plt

pplt = plot_wave_real(data, meta, meta['Nt'],
                      n = 10, cmapname='plasma').show()




# %% ======================================================
#  PLOT |  Real & complex part of wavefunction
# =========================================================

def plot_wave_both(data, meta, tmaxidx, n,
              basp = 40,
              figsize = (12, 6), cmapname='plasma'):
    cmap = plt.get_cmap(cmapname)
    fig, ax = plt.subplots(1, 2, figsize=figsize)

    timevals = []
    for i in range(0, tmaxidx, int(tmaxidx/n) ):
        ax[0].plot(np.linspace(-meta['xmax'], meta['xmax'], num=meta['Nx']),
                 np.real(data[i,:]),
                 color = cmap(i/(tmaxidx-1))
                )
        ax[1].plot(np.linspace(-meta['xmax'], meta['xmax'], num=meta['Nx']),
                 np.imag(data[i,:]),
                 color = cmap(i/(tmaxidx-1))
                )
        timevals.append(meta['t'][i])

    norm = mpl.colors.Normalize(vmin=timevals[0], vmax=timevals[-1])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    ax[0].set_title("$Re(\psi)$")
    ax[1].set_title("$Im(\psi)$")
    ax[0].set_xlabel("$x$", labelpad=2)
    ax[1].set_xlabel("$x$", labelpad=2)
    #plt.xlabel("$x$", labelpad=2)
    #plt.ylabel("$Re(\psi)$")
    plt.colorbar(sm, ticks=[round(x) for x in timevals], ax=ax,
                 orientation='horizontal',
                 aspect=basp,
                 label='$t$')
    return plt


pplt = plot_wave_both(data, meta, 5010,
                      n = 10,
                      cmapname='plasma')

if do_writeplot: 
    pplt.savefig( fname_plots.format("optimal_wave"), transparent=True)



# %% ======================================================
#  PLOT |  fit max 
# =========================================================

block_figsize, block_aspect = (12,6), 0.6
block_lim = (-2,6)

pplt = plot_density_timeline(data, meta, cmap='viridis',
                             ylim = block_lim,
                             figsize = block_figsize, aspect = block_aspect)


# %% ======================================================
#  DATA | load other data to show fit result, tau = 4
# =========================================================

data, meta = read_fortran_hw5('data/psi_10001-10000-4.0.dat')

pplt = plot_density_timeline(data, meta, cmap='viridis',
                             ylim = block_lim,
                             figsize = block_figsize, aspect = block_aspect)
add_fit(data, meta)
plt.xticks([]); plt.xlabel(""); # customization for pdf

if do_writeplot: 
    pplt.savefig( fname_plots_png.format("psi_tau-4_fit"), transparent=True)



# %% ======================================================
#  DATA | load other data to show fit result, tau = 8
# =========================================================

data, meta = read_fortran_hw5('data/psi_10001-10000-8.0.dat')

pplt = plot_density_timeline(data, meta, cmap='viridis',
                             ylim = block_lim,
                             figsize = block_figsize, aspect = block_aspect)
add_fit(data, meta)
plt.xticks([]); plt.xlabel(""); # customization for pdf

if do_writeplot: 
    pplt.savefig( fname_plots_png.format("psi_tau-8_fit"), transparent=True)


# %% ======================================================
#  DATA | load other data to show fit result, tau = 15
# =========================================================

data, meta = read_fortran_hw5('data/psi_10001-10000-15.0.dat')

pplt = plot_density_timeline(data, meta, cmap='viridis',
                             ylim = block_lim,
                             figsize = block_figsize, aspect = block_aspect)
add_fit(data, meta)

if do_writeplot: 
    pplt.savefig( fname_plots_png.format("psi_tau-15_fit"), transparent=True)












# %% ======================================================
#  DATA | load subsampled dataset
# =========================================================

data, meta = read_fortran_hw5('data/psi_10001-200-10.0.dat')
print_metadata(meta)




# %%

pplt = plot_density_timeline(data, meta)
if do_writeplot: 
    pplt.savefig( fname_plots_png.format("psi_subsampled_density"), transparent=True)

# %%

pplt = plot_wave_both(data, meta, 200, n = 3,
                      cmapname='plasma')
if do_writeplot:
    pplt.savefig( fname_plots.format("psi_subsampled_wave"), transparent=True)






# %% ======================================================
#  DATA | load psi_3 test data
# =========================================================

data, meta = read_fortran_hw5('data/psi3_10001-10000-10.0.dat')
print_metadata(meta)


# %%

pplt = plot_density_timeline(data, meta, ylim=(-4,6), aspect=0.8)
if do_writeplot: 
    pplt.savefig( fname_plots_png.format("psi3_density"), transparent=True)

# %%

pplt = plot_wave_both(data, meta, 3000, n = 3,
                      cmapname='plasma')
if do_writeplot: 
    pplt.savefig( fname_plots.format("psi3_wave"), transparent=True)







# %% ======================================================
#  DATA | load psi_10 test data
# =========================================================

data, meta = read_fortran_hw5('data/psi10_10001-10000-10.0.dat')
print_metadata(meta)


# %%

pplt = plot_density_timeline(data, meta, ylim=(-6,8), aspect=0.57)
if do_writeplot: 
    pplt.savefig( fname_plots_png.format("psi10_density"), transparent=True)

# %%

pplt = plot_wave_both(data, meta, 3000, n = 3,
                      cmapname='plasma')
if do_writeplot: 
    pplt.savefig( fname_plots.format("psi10_wave"), transparent=True)





# %% ======================================================
#  DATA | load degraded dataset (with problems at margin)
# =========================================================

data, meta = read_fortran_hw5('data/psi_10001-10000-1.0.dat')
print_metadata(meta)


# %%

pplt = plot_density_timeline(data, meta, ylim=(-3,10), aspect=0.8)
if do_writeplot: 
    pplt.savefig( fname_plots_png.format("psi_tau1_density"), transparent=True)

# %%

pplt = plot_wave_both(data, meta, 10000, n = 5,
                      cmapname='plasma')
if do_writeplot: 
    pplt.savefig( fname_plots.format("psi_tau1_wave"), transparent=True)





# %% ======================================================
#  DATA | load omega=3.0 dataset (to show oscillation component)
# =========================================================


data, meta = read_fortran_hw5('data/psi_omega3_10001-10000-10.0.dat')
print_metadata(meta)

pplt = plot_density_timeline(data, meta, cmap='viridis',
                             ylim = block_lim,
                             figsize = block_figsize, aspect = block_aspect)
if do_writeplot: 
    pplt.savefig( fname_plots_png.format("omega3"), transparent=True)

pplt = plot_linearized_residual(data, meta)

if do_writeplot: 
    pplt.savefig( fname_plots.format("omega3_res"), transparent=True)

# %%
