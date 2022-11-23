# %% general import ---------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

do_writeplot = True
do_forceclean = True
fname_plots = "img/{}.pdf"

# %% set the plotting env ---------------------------------

useLaTeX = True  # if True, use LaTeX backend
if useLaTeX:
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)

plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (7,4)

# %% import data ------------------------------------------
files = {
        'diag' : "data/exr3_diag.txt",
        'herm' : "data/exr3_herm.txt"
}

data = {}
for k, v in files.items():
    print("loading file {} ({})".format(v, k) )
    spacing = np.loadtxt(v)
    print( " imported {} values".format(len(spacing)) )
    data[k] = spacing
# %% other defs

# tools ---------------------------------------------------
def make_bins(x, step, n):
    hist, bins = np.histogram(x, bins = [ step*k for k in range(0,n)], density=True)
    bins = bins[:-1] + step/2 # take the middle point for each bin, for fitting later
    return hist, bins

# fitting functions ---------------------------------------
def generic_fit(x, a, b, alpha, beta):
    return a * (x**alpha) * np.exp(-b *(x**beta))

def wigner_surmise(x, a, b):
    return a * (x**2) * np.exp(-b *(x**2))

def quasi_poisson(x, a, b, beta):
    return a * np.exp(-b *(x**beta))

def poisson(x, a, b):
    return a * np.exp(-b*x)


# %%  setting default parameters --------------------------
plt_xrange_general = (-0.07, 4.07) # set global xlim of plots
step = 0.1   # default bin width value
maxn = 200   # max bin product factor


# %% ======================================================
#  PLOT 0 |  general histogram
# =========================================================

# binning & plot
hist, bins = make_bins(data['herm'], step, maxn)
plt.bar(bins, hist, width=step, color='goldenrod', edgecolor='k', linewidth=1, label='Hermitian', alpha = 0.7)
hist, bins = make_bins(data['diag'], step, maxn)
plt.bar(bins, hist, width=step, color='mediumorchid', edgecolor='k', linewidth=1, label="Diagonal", alpha = 0.7)

# set plot properties
plt.xlim( plt_xrange_general )
plt.xlabel("$s$");  plt.ylabel("$P(s)$"); 
plt.locator_params(axis='y', nbins=6)
plt.legend()
plt.tight_layout()
if do_writeplot: plt.savefig( fname_plots.format("exr3-general"), transparent=True)
if do_forceclean: plt.clf()




# %% ======================================================
#  PLOT 1 |  Hermitian matrices, full fit
# =========================================================

plt_xrange_general = (-0.07, 3.07)

# binning
hist, bins = make_bins(data['herm'], step, maxn)
# fitting
popt, pcov = curve_fit(generic_fit, bins, hist)
print('fit coeffs :', popt)
print('  with err :', np.sqrt(np.diag(pcov)) )

# plot fit
xx =  np.arange(0.001,plt_xrange_general[1],0.001)
plt.plot(xx, generic_fit(xx, *popt), 'r', linewidth=3, label="full fit")
# plot bars
plt.bar(bins,hist, width=step, color='moccasin', edgecolor='k', linewidth=1, label='Hermitian' )

# set plot properties
plt.xlim( plt_xrange_general )
plt.xlabel("$s$");  plt.ylabel("$P(s)$"); 
plt.locator_params(axis='y', nbins=6)
plt.legend()
plt.tight_layout()
if do_writeplot: plt.savefig( fname_plots.format("herm_full-fit"), transparent=True)
if do_forceclean: plt.clf()



# %% ======================================================
#  PLOT 2 |  Diagonal matrices, full fit
# =========================================================

# binning
hist, bins = make_bins(data['diag'], step, maxn)
# fitting
popt, pcov = curve_fit(generic_fit, bins, hist)
print('fit coeffs :', popt)
print('  with err :', np.sqrt(np.diag(pcov)) )

# plot fit
xx =  np.arange(0.001,plt_xrange_general[1],0.001)
plt.plot(xx, generic_fit(xx, *popt), 'r', linewidth=3, label='full fit')
# plot bars
plt.bar(bins,hist, width=step, color='thistle', edgecolor='k', linewidth=1, label='Diagonal')

# set plot properties
plt.xlim( plt_xrange_general )
plt.xlabel("$s$");  plt.ylabel("$P(s)$"); 
plt.locator_params(axis='y', nbins=6)
plt.legend()
plt.tight_layout()
if do_writeplot: plt.savefig( fname_plots.format("diag_full-fit"), transparent=True)
if do_forceclean: plt.clf()



# %% ======================================================
#  PLOT 3 |  Hermitian matrices, coefficients only fit
# =========================================================

step = 0.05   # default bin width value
maxn = 400    # max bin product factor

# binning
hist, bins = make_bins(data['herm'], step, maxn)
# fitting
popt, pcov = curve_fit(wigner_surmise, bins, hist)
print('fit coeffs :', popt)
print('  with err :', np.sqrt(np.diag(pcov)) )

# plot fit
xx =  np.arange(0.001, plt_xrange_general[1], 0.001)
plt.plot(xx, wigner_surmise(xx, *popt), 'r', linewidth=3, label='coeff fit')
plt.plot(xx, wigner_surmise(xx, 32/(np.pi**2),4/np.pi), 'dodgerblue', linewidth=2, label='Wigner surmise')
# plot bars
plt.bar(bins,hist, width=step, color='moccasin', edgecolor='k', linewidth=1, label='Hermitian')

# set plot properties
plt.xlim( plt_xrange_general )
plt.xlabel("$s$");  plt.ylabel("$P(s)$"); 
plt.locator_params(axis='y', nbins=6)
plt.legend()
plt.tight_layout()
if do_writeplot: plt.savefig( fname_plots.format("herm_coeff-fit"), transparent=True)
if do_forceclean: plt.clf()



# %% ======================================================
#  PLOT 4 |  Diagonal matrices, coefficients only fit
# =========================================================

# binning
hist, bins = make_bins(data['diag'], step, maxn)
# fitting
popt, pcov = curve_fit(poisson, bins, hist)
print('fit coeffs :', popt)
print('  with err :', np.sqrt(np.diag(pcov)) )

# plot fit
xx =  np.arange(0.001, plt_xrange_general[1], 0.001)
plt.plot(xx, poisson(xx, *popt), 'r', linewidth=3, label='poisson fit')
# plot bars
plt.bar(bins,hist, width=step, color='thistle', edgecolor='k', linewidth=1, label='Diagonal')

# set plot properties
plt.xlim( plt_xrange_general )
plt.xlabel("$s$");  plt.ylabel("$P(s)$"); 
plt.locator_params(axis='y', nbins=6)
plt.legend()
plt.tight_layout()
if do_writeplot: plt.savefig( fname_plots.format("diag_coeff-fit"), transparent=True)
if do_forceclean: plt.clf()

