# ========================================
#  QIC  >>  assignment 4
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 23 November 2022
# ========================================
#
# This script accounts for the analysis of
# homework 4 data.
#
#  Requires :  data/*.eigvec  (binary)
#              data/*.eigval  (txt)
#               (in pairs)
#
# ========================================



# %% general import ---------------------------------------

import numpy as np
import numpy.polynomial.hermite as Herm
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import pandas as pd
import os
from scipy.io import FortranFile


# [!]  you might want to test  'tri/'  and  'penta/'
SUBSET = 'tri/' 

do_writeplot = True
do_forceclean = False
img_prefix = "tex/img/" + SUBSET
fname_plots = img_prefix + "{}.pdf"


# %% set the plotting env ---------------------------------

useLaTeX = True  # if True, use LaTeX backend
if useLaTeX:
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)

plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (9,6)



# %% ======================================================
#  IMPORT DATA
# =========================================================

# target data directory
PREFIX = 'data/' + SUBSET

# target file extensions
POSTFIX = { 'val': '.eigval',   # extension of eigenvalue file
            'vec' : '.eigvec'}  # extension of eigenvector file

# look for files in target directory
data_files = [f.strip(POSTFIX['val']) for f in os.listdir(PREFIX) if f.endswith(POSTFIX['val'])]
print( "found {} records".format(len(data_files)) )

# to read Fortran binary files
def read_real8_matrix(fname : str) -> np.matrix:
    """Read binary real8 matrix file"""
    f = FortranFile( fname, 'r')
    shape = f.read_ints(np.int32)
    mx = f.read_reals(float)
    return mx.reshape(shape, order="F")

# load data from target directory
temp_df = []
for dfile in sorted(data_files):
    # parse string parameters
    pp = dfile.split('-')
    N = int( pp[0] );   xmax = float(pp[1]);   omega = float(pp[2]); 
    
    # set prefix
    dfile = PREFIX + dfile              
    print( "loading {}".format(dfile) )

    # load data
    eigval = np.loadtxt(dfile + POSTFIX['val'])
    eigvec = read_real8_matrix(dfile + POSTFIX['vec']).T # to reaad from binary
    # np.genfromtxt(dfile + POSTFIX['vec'],  delimiter=",", comments='#') # to read from csv

    # add to temp dataframe + preprocessing
    tdf = pd.DataFrame(data = { 'N': [N], 'xmax': [xmax], 'omega': [omega] },
                       columns=['N', 'xmax', 'omega', 'dx', 'eigval', 'eigvec'])
    tdf['eigval'] = pd.Series(dtype="object")
    tdf['eigvec'] = pd.Series(dtype="object")

    tdf['dx'] = [ 2*xmax/(N-1) ]
    tdf.at[0,'eigval'] = eigval
    tdf.at[0,'eigvec'] = eigvec
    temp_df.append(tdf)

df = pd.concat(temp_df)
df = df.sort_values(by=['N'])
df




# %% ======================================================
#  utilities
# =========================================================

# Remark: at this link
#   https://chem.libretexts.org/Ancillary_Materials/Interactive_Applications/Jupyter_Notebooks/Quantum_Harmonic_Oscillators_-_Plotting_Eigenstates_(Python_Notebook)
# there is some cool stuff about harmonic oscillators

class harmonic_theory():
    # for hbar = 1, m = 1

    def __init__(self, omega):
        self.omega = omega

    def hermite(self, x, n):
        xi = np.sqrt(self.omega)*x
        herm_coeffs = np.zeros(n+1)
        herm_coeffs[n] = 1
        return Herm.hermval(xi, herm_coeffs)

    def eigenstates(self, x, n : int):
        coeff = (self.omega/np.pi)**(0.25)/np.sqrt( (2.**n) * np.math.factorial(n))
        wavef = coeff * np.exp(- self.omega*(x**2) / 2) * self.hermite(x, n)
        return wavef

    def eigenvalues(self, n):
        return self.omega*(n + .5)


def get_wavef(dsel: dict, n : int):
    """This function retrieves an eigenstate from the dataset."""

    # unwrap parameters
    xmax = dsel['xmax']

    x = np.linspace(-xmax, xmax, num=dsel['N'])
    # [!] deprecated due to bug in accumulation precision ---
    #dx = 2*xmax/(dsel['N']-1)
    #x = np.arange(-xmax, xmax+dx, dx ) 
    # -------------------------------------------------------

    # normalize the wavefunction
    wavef = dsel['eigvec'][n,:]
    wavef = wavef/np.linalg.norm(wavef, ord=2)  # just to be sure, normalize

    # remark: the wavefunction is automatically flipped to
    #         be mostly positive for x > 0
    return x, wavef*np.sign( np.sum(wavef[int(len(wavef)/2):] ) )

def error_eigval(dsel : dict):
    """Compute relative error of eigenvalues wrt theoretically expected values.

    Parameters
    ----------
    dsel : dict
        The dictionary corresponding to a dataset row.
    
    Returns
    -------
    array
        The relative error of each eigenvalue associated to psi_n wrt theory
    """

    this_eigval = dsel['eigval']

    ht = harmonic_theory(dsel['omega'])
    theo_eigval = ht.eigenvalues( n = np.arange( len(this_eigval) ) )

    return np.abs(this_eigval - theo_eigval)/theo_eigval





# %% ======================================================
#  PLOT |  Eigenstates density (psi^2)
# =========================================================

# select a reference record
selection = df[df['N'] == 1001].to_dict('records')[0]
xmax = selection['xmax']
# max n-1 of wavefunction to plot
max_ord = 21
# visual amplification factor for this plot
ampl = 30 
plt.figure(figsize=(7,9))

for nn in range(max_ord):
    x, wf = get_wavef(selection, n = nn)
    plt.plot(x, ampl*(wf**2)+nn, color='k')
    plt.fill_between(x, ampl*(wf**2)+nn, nn,
                     color='k', alpha=0.3)

plt.ylim([-0.5, max_ord-0.2])
plt.xlim([-10.0,10.0])
plt.yticks(ticks=[n for n in range(max_ord)])
plt.title(r"$|\psi_k|^2$", y=1.02)
plt.xlabel(r"$x$")
plt.ylabel(r"$k$", rotation=0, labelpad = 15)
plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigstates_density"), transparent=True)
if do_forceclean: plt.clf()




# %% ======================================================
#  PLOT |  Eigenstates density (psi^2) with colormap
# =========================================================

# select a reference record
selection = df[df['N'] == 1001].to_dict('records')[0]

# max n-1 of wavefunction to plot
max_ord = 15

# create matrix for heatmap
collection = []
for nn in range(max_ord):
    x, wf = get_wavef(selection, n = nn)
    collection.append(wf**2)

plt.figure(figsize=(7,9))
plt.imshow(collection, cmap='hot',
           interpolation='none', aspect=1.7,
           extent=[-selection['xmax'],selection['xmax'], max_ord-.5,-.5])

plt.ylim([-0.5, max_ord-0.5])
plt.xlim([-10.0, 10.0])
plt.yticks(ticks=[k for k in range(max_ord)])
plt.title(r"$|\psi_k|^2$", y=1.02)
plt.xlabel(r"$x$")
plt.ylabel(r"$k$", rotation=0, labelpad = 15)
plt.tight_layout()

if do_writeplot: 
    plt.savefig( img_prefix + "eigstates_density_cmap.png", transparent=True)
if do_forceclean: plt.clf()



# %% ======================================================
#  PLOT |  Eigenstates (psi_n)
# =========================================================

# select a reference record
selection = df[df['N'] == 1001].to_dict('records')[0]
# max n-1 of wavefunction to plot
max_ord = 11

plt.figure(figsize=(8,8))

ht = harmonic_theory(selection['omega'])
ampl = 3 # visual amplification factor for this plot

for nn in range(max_ord):
    x, wf = get_wavef(selection, n = nn)
    
    plt.plot(x, ampl*wf+nn, color='k')
    plt.fill_between(x, ampl*wf+nn, nn,
                     color='k', alpha=0.3)

plt.ylim([-0.5, max_ord-0.2])
plt.xlim([-10.0,10.0])
plt.yticks(ticks=[n for n in range(max_ord)])
plt.title(r"$\psi_k$", y=1.01)
plt.xlabel(r"$x$")
plt.ylabel(r"$k$", rotation=0, labelpad = 15)
plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigstates_raw"), transparent=True)
if do_forceclean: plt.clf()




# %% ======================================================
#  PLOT |  Eigenvalues relative error (function of order)
# =========================================================


plt.figure(figsize=(9,8))

N_select = [201,401,801,1001,2001,4001,8001,10001]
max_range = min(N_select)
colors = pl.cm.plasma(np.linspace(0.85,0,len(N_select)))
for idx, NN in enumerate(N_select):
    ss = df[df['N'] == NN].to_dict('records')[0]
    plt.plot( error_eigval(ss)[:max_range],
              #linestyle='--', marker='x',
              linewidth=2.5, #markersize=2,
              label="${}$".format(NN), color=colors[idx] )

# add 1% line
plt.axhline(y = 0.01, color = 'g', linestyle = '--', linewidth=1)
plt.text(182, 0.011, r"$1\%\;\textrm{error}$", size=16)

leg = plt.legend(title=r"$N$", prop={'size': 18}, ncol=4)
for legobj in leg.legendHandles:
    legobj.set_linewidth(4)

plt.ylabel(r"$\varepsilon_{rel}^{(k)}$",rotation=0, y=0.92, labelpad=16)
plt.xlabel(r"$k$")
plt.yscale("log")
plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigval_relative_err"), transparent=True)
if do_forceclean: plt.clf()



# %% ======================================================
#  PLOT |  Eigenvalues (raw) (function of order)
# =========================================================

max_range = 200 - 1
N_select = [1001,2001,4001,8001]
colors = pl.cm.plasma(np.linspace(0.85,0,len(N_select)))
for idx, NN in enumerate(N_select):
    ss = df[df['N'] == NN].to_dict('records')[0]
    plt.plot( np.arange(1,max_range+1), ss['eigval'][1:max_range+1],
              #linestyle='--', marker='x',
              #linewidth=1, #markersize=2,
              label="${}$".format(NN), color=colors[idx] )

ht = harmonic_theory(ss['omega'])
theo_eigval = ht.eigenvalues( n = np.arange(1,max_range+1) )
plt.plot( np.arange(1,max_range+1), theo_eigval, label=r"$\textrm{theory}$",
          linestyle='--', color='r')

plt.legend(title=r"$N$", prop={'size': 20})
plt.ylabel(r"$\varepsilon_{k}$",rotation=0, labelpad=14)
plt.xlabel(r"$k$", horizontalalignment='right', x=1.0)
plt.yscale("log")
plt.xscale("log")
plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigval_raw"), transparent=True)
if do_forceclean: plt.clf()



# %% ======================================================
#  PLOT |  Eigenvalues relative err (function of N)
#       |  for first 4 eigenvalues
# =========================================================
n_eigenvs = 4

sel_eigvals = [k for k in range(n_eigenvs)]
errors=[]
Ns = []
dxs = []
for ss in df.to_dict('records'):
    errors.append( error_eigval(ss)[sel_eigvals] )
    Ns.append(ss['N'])
    dxs.append(ss['dx'])

fig, ax = plt.subplots(constrained_layout=True)
colors = pl.cm.plasma(np.linspace(0,0.85,n_eigenvs))
for ee, vv in zip(sel_eigvals, np.matrix(errors).T):
    ax.plot(Ns, vv.T, label="${}$".format(ee), color=colors[ee])

ax.set_ylabel(r"$\varepsilon_{rel}^{(k)}$", rotation=0, labelpad=28)
ax.set_xlabel(r"$N$")
ax.set_yscale("log")


def forward(x): # from primary to secondary axis
    return 2*xmax/(x-1)

def inverse(x):
    return (2*xmax/x)+1

xtick = ax.get_xticks()[2:-1]
print('ref xticks for numbers in bottom axis', forward(xtick) )

secax = ax.secondary_xaxis('top', functions=(forward, inverse))
secax.set_xlabel(r"$dx$")
secax.set_xticks( [0.1,0.01,0.005,0.003,0.002] ) # to be manually tweaked!

plt.legend(title="$k$", loc='upper right', ncols=n_eigenvs)
plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigval_relative_err_N"), transparent=True)
if do_forceclean: plt.clf()




# %% ======================================================
#  PLOT |  Eigenfunctions plot
# =========================================================

def compare_eigenfunctions(dsel : dict, n : int):
    x, wf = get_wavef(dsel, n = n)

    ht = harmonic_theory(dsel['omega'])
    wftheo = ht.eigenstates(x, n)
    wftheo = wftheo/np.linalg.norm(wftheo, ord=2)

    plt.plot(x, wf, color='k', label='numerical')
    plt.plot(x, wftheo, color='r', linestyle='--', label='theory')
    plt.legend()
    plt.ylabel(r"$\psi_{%d}$" % (n),rotation=0, labelpad=8)
    plt.xlabel(r"$x$")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.tight_layout()
    return plt

# plot psi_1
pplt = compare_eigenfunctions(selection, 1)
if do_writeplot:
    pplt.savefig( fname_plots.format("eigstates_ex01"), transparent=True)
    pplt.clf()

# plot psi_15
pplt = compare_eigenfunctions(selection, 15)
if do_writeplot: 
    pplt.savefig( fname_plots.format("eigstates_ex15"), transparent=True)
    pplt.clf()

# plot psi_40
pplt = compare_eigenfunctions(selection, 40)
if do_writeplot: 
    pplt.savefig( fname_plots.format("eigstates_ex40"), transparent=True)
    pplt.clf()




# %% ======================================================
#  PLOT |  Eigenfunction comulative error (function of order)
# =========================================================

def err_eigenfunctions(dsel : dict, n : int):
    x, wf = get_wavef(dsel, n = n)

    ht = harmonic_theory(dsel['omega'])
    wftheo = ht.eigenstates(x, n)
    wftheo = wftheo/np.linalg.norm(wftheo, ord=2)

    return np.sum( np.abs(wftheo - wf) )/dsel['N']


max_range = 100
N_select = [401,801,1001,2001,4001,8001,10001]
colors = pl.cm.plasma(np.linspace(0.85,0,len(N_select)))
for idx, NN in enumerate(N_select):
    ss = df[df['N'] == NN].to_dict('records')[0]
    yplot = [ err_eigenfunctions(ss, i) for i in range(max_range) ]

    plt.plot( yplot,
              linestyle='--', marker='.',
              linewidth=1, #markersize=5,
              label="${}$".format(NN), color=colors[idx] )

leg = plt.legend(title=r"$N$", prop={'size': 20}, ncol=3)
for legobj in leg.legendHandles:
    legobj.set_markersize(22)

plt.ylabel(r"$\delta\psi_k$",rotation=0, labelpad=23)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlabel(r"$k$")
plt.yscale("log")
plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigstates_err_comulative"), transparent=True)
if do_forceclean: plt.clf()




# %% exit

print('End of script')