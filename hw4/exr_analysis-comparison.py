# ========================================
#  QIC  >>  assignment 4
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 25 November 2022
# ========================================
#
# This script accounts for the analysis of
# homework 4 data. In particular, here
# I wish to compare the tridiagonal FDM
# with the pentadiagonal.
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
from matplotlib.legend_handler import HandlerTuple
import pandas as pd
import os
from scipy.io import FortranFile



do_writeplot = True
do_forceclean = False
img_prefix = "tex/img/compare/"
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
PREFIX = 'data/'

# target file extensions
POSTFIX = { 'val': '.eigval',   # extension of eigenvalue file
            'vec' : '.eigvec'}  # extension of eigenvector file

# to read Fortran binary files
def read_real8_matrix(fname : str) -> np.matrix:
    """Read binary real8 matrix file"""
    f = FortranFile( fname, 'r')
    shape = f.read_ints(np.int32)
    mx = f.read_reals(float)
    return mx.reshape(shape, order="F")


def load_dataset(_subpath):
    # look for files in target directory
    data_files = [f.strip(POSTFIX['val']) for f in os.listdir(PREFIX + _subpath) if f.endswith(POSTFIX['val'])]
    print( "found {} records".format(len(data_files)) )

    # load data from target directory
    temp_df = []
    for dfile in sorted(data_files):
        # parse string parameters
        pp = dfile.split('-')
        N = int( pp[0] );   xmax = float(pp[1]);   omega = float(pp[2]); 

        # set prefix
        dfile = PREFIX + _subpath + dfile              
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

    return pd.concat(temp_df)

df_tri = load_dataset('tri/')
df_penta = load_dataset('penta/')




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


class HandlerTupleVertical(HandlerTuple):
    """Plots all the given Lines vertical stacked."""

    def __init__(self, **kwargs):
        """Run Base Handler."""
        HandlerTuple.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        """Create artists (the symbol) for legend entry."""
        # How many lines are there.
        numlines = len(orig_handle)
        handler_map = legend.get_legend_handler_map()

        # divide the vertical space where the lines will go
        # into equal parts based on the number of lines
        height_y = (height / numlines)

        leglines = []
        for i, handle in enumerate(orig_handle):
            handler = legend.get_legend_handler(handler_map, handle)

            legline = handler.create_artists(legend, handle,
                                             xdescent,
                                             (2*i + 1)*height_y,
                                             width,
                                             2*height,
                                             fontsize, trans)
            leglines.extend(legline)

        return leglines


# %% ======================================================
#  requirements
# =========================================================

selection = df_tri[df_tri['N'] == 1001].to_dict('records')[0]
xmax = selection['xmax'] # Hp: homogeneous xmax



# %% ======================================================
#  PLOT |  Eigenvalues relative error (function of order)
# =========================================================


fig, ax = plt.subplots(figsize=(9,8))

N_select = [201,601, 1001,2001,4001,10001]
max_range = min(N_select)

colors = pl.cm.winter(np.linspace(1,0,len(N_select)))
for idx, NN in enumerate(N_select):
    ss = df_tri[df_tri['N'] == NN].to_dict('records')[0]
    ax.plot( error_eigval(ss)[:max_range],
              linestyle='--', #marker='x',
              linewidth=2.5, #markersize=2,
              label="${}$".format(NN), color=colors[idx] )


max_range = min(N_select)
colors = pl.cm.autumn(np.linspace(0.80,0,len(N_select)))
for idx, NN in enumerate(N_select):
    ss = df_penta[df_penta['N'] == NN].to_dict('records')[0]
    ax.plot( error_eigval(ss)[:max_range],
              #linestyle='--', marker='x',
              linewidth=2.5, #markersize=2,
              label="${}$".format(NN), color=colors[idx] )

handles, labels = ax.get_legend_handles_labels()

# unwrapping the legend
triHandles = handles[:len(N_select)]
triArtist = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
pentaHandles = handles[len(N_select):]
pentaArtist = plt.Line2D((0,1),(0,0), color='k')
labels = labels[:len(N_select)]

legg = ax.legend( [ tuple(x) for x in zip(triHandles, pentaHandles)],
                  labels , title="$N$", ncol = 3,
           handler_map={tuple: HandlerTupleVertical(ndivide=None)},
           loc='lower center'
)

ax.legend( [triArtist, pentaArtist], ['tridiagonal', 'pentadiagonal'],
           loc='upper center', bbox_to_anchor=(0.5, 1.10),
           ncol = 2, fontsize='small',
           #title=r"finite difference method"
)

ax.add_artist(legg)

plt.ylabel(r"$\varepsilon_{rel}^{(k)}$",rotation=0, labelpad=16)
plt.xlabel(r"$k$")
plt.yscale("log")
plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigval_relative_err"), transparent=True)
if do_forceclean: plt.clf()



# %% ======================================================
#  PLOT |  Eigenvalues relative err (function of N)
#       |  for first 4 eigenvalues
# =========================================================
n_eigenvs = 4
sel_eigvals = [k for k in range(n_eigenvs)]
fig, ax = plt.subplots(constrained_layout=True)
colors = pl.cm.plasma(np.linspace(0.80,0,n_eigenvs))

errors=[]
Ns = []
dxs = []
for ss in df_penta.sort_values(by=['N']).to_dict('records'):
    errors.append( error_eigval(ss)[sel_eigvals] )
    Ns.append(ss['N'])
    dxs.append(ss['dx'])

for ee, vv in zip(sel_eigvals, np.matrix(errors).T):
    ax.plot(Ns, vv.T, label="${}$".format(ee), color=colors[ee])

errors=[]
Ns = []
dxs = []
for ss in df_tri.sort_values(by=['N']).to_dict('records'):
    errors.append( error_eigval(ss)[sel_eigvals] )
    Ns.append(ss['N'])
    dxs.append(ss['dx'])

for ee, vv in zip(sel_eigvals, np.matrix(errors).T):
    ax.plot(Ns, vv.T, label="${}$".format(ee), color=colors[ee], linestyle='--')

handles, labels = ax.get_legend_handles_labels()

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


# unwrapping the legend
triArtist = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
pentaArtist = plt.Line2D((0,1),(0,0), color='k')

legg = ax.legend( handles[:-4], labels[:-4] , title="$N$", ncol = n_eigenvs,
           handler_map={tuple: HandlerTupleVertical(ndivide=None)},
           loc='upper right',  prop={'size': 17}
)

ax.legend( [triArtist, pentaArtist], ['tridiagonal', 'pentadiagonal'],
           loc='lower left',
           ncol = 1, fontsize='small', prop={'size': 17},
)

ax.add_artist(legg)

plt.tight_layout()

if do_writeplot: 
    plt.savefig( fname_plots.format("eigval_relative_err_N"), transparent=True)
if do_forceclean: plt.clf()





# %% ======================================================
#  PLOT |  Eigenfunction comulative error (function of order)
# =========================================================

def err_eigenfunctions(dsel : dict, n : int):
    x, wf = get_wavef(dsel, n = n)

    ht = harmonic_theory(dsel['omega'])
    wftheo = ht.eigenstates(x, n)
    wftheo = wftheo/np.linalg.norm(wftheo, ord=2)

    return np.sum( np.abs(wftheo - wf) )/dsel['N']

fig, ax = plt.subplots(figsize=(9,8))

max_range = 100
N_select = [401,801,1001,2001,6001,10001]

colors = pl.cm.winter(np.linspace(1,0,len(N_select)))
for idx, NN in enumerate(N_select):
    ss = df_tri[df_tri['N'] == NN].to_dict('records')[0]
    yplot = [ err_eigenfunctions(ss, i) for i in range(max_range) ]

    ax.plot( yplot, linestyle='--', linewidth=2.5,
              label="${}$".format(NN), color=colors[idx] )


colors = pl.cm.autumn(np.linspace(0.80,0,len(N_select)))
for idx, NN in enumerate(N_select):
    ss = df_penta[df_penta['N'] == NN].to_dict('records')[0]
    yplot = [ err_eigenfunctions(ss, i) for i in range(max_range) ]

    ax.plot( yplot, linewidth=2.5,
              label="${}$".format(NN), color=colors[idx] )


# unwrapping the legend
handles, labels = ax.get_legend_handles_labels()

triHandles = handles[:len(N_select)]
triArtist = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
pentaHandles = handles[len(N_select):]
pentaArtist = plt.Line2D((0,1),(0,0), color='k')
labels = labels[:len(N_select)]

legg = ax.legend( [ tuple(x) for x in zip(triHandles, pentaHandles)],
                  labels , title="$N$", ncol = 3,
           handler_map={tuple: HandlerTupleVertical(ndivide=None)},
           loc='lower right'
)

ax.legend( [triArtist, pentaArtist], ['tridiagonal', 'pentadiagonal'],
           loc='upper center', bbox_to_anchor=(0.5, 1.10),
           ncol = 2, fontsize='small',
           #title=r"finite difference method"
)

ax.add_artist(legg)


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
# %%
