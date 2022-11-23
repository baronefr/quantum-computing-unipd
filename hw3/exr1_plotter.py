#!/bin/python3

# ========================================
#  QIC  >>  assignment 3 / exr 1
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 09 November 2022
# ========================================


import matplotlib.pyplot as plt
import pandas as pd
from scipy import odr
import numpy as np

# import from other module the configurations
from exr1 import *
print("looking for opti    : ", OPTI_LEVELS)
print("looking for dataset : ", OPTI_BENCHFILES)


# I/O management ---------------------------------------------------
PREFIX_OUTPUT = "./img/"
fname_plots      = PREFIX_OUTPUT + "{}" + '.pdf'
do_writeplot = True   # if True, plots are written in file
# ------------------------------------------------------------------


# set the plotting env ---------------------------------------------

useLaTeX = True  # if True, use LaTeX
if useLaTeX:
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)

plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (10,6)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# import dataset as dataframe
# ------------------------------------------------------------------
temp_df = []
for bfile, op in zip(OPTI_BENCHFILES, OPTI_LEVELS):
    print("retrieve ", bfile)
    tdf = pd.read_csv( bfile,
                       names = ['dim','loop RbC','loop CbR','loop JKI','matmul'] )
    tdf['opti'] = int(op)
    temp_df.append(tdf)

df = pd.concat(temp_df)
df


# make plot dir if needed
if do_writeplot:
    if not os.path.exists(PREFIX_OUTPUT):
        print('creating output directory')
        os.makedirs(PREFIX_OUTPUT)


# ------------------------------------------------------------------
# PLOT 1: fit CPU time with order 3 polynomial
# ------------------------------------------------------------------
print("making plot 1")
no_opti = df[ df['opti']==0 ]

ax = no_opti.plot( x = 'dim', y = ['loop RbC','loop CbR','matmul'],  
                label = ['RbC','CbR','matmul'], #figsize=(7,6), 
                linewidth=3 )
ax.set( xlabel='matrix dimension $N$', ylabel='time $[s]$',
        title = 'gfortran: no compiler optimization')

# fit O(n^3) ------------
output = np.polyfit(no_opti['dim'], no_opti['loop RbC'], 3)
poly = np.poly1d(output)
print(' fit poly3 coeffs: ', output)

x = np.arange(0, 1.05*no_opti['dim'].max(), 10)
plt.plot(x, poly(x), ':', label="$O(N^3)$ fit", linewidth=2, zorder=3)
plt.legend()
plt.tight_layout()

if do_writeplot:
    fname = fname_plots.format("poly3_fit")
    print("writing plot", fname)
    plt.savefig( fname, transparent=True)




# ------------------------------------------------------------------
# PLOT 2: compare the CPU time with different optimizations
# ------------------------------------------------------------------
print("making plot 2")

# reference: matmul & worse timing for any loop
ax = df[ df['opti'] == 5 ].plot(x='dim', y = 'matmul', label='matmul', color='k')
ax = df[ df['opti'] == 0 ].plot(x='dim', y='loop RbC', label='RbC -O0', color='red', ax=ax, linewidth=2) # these curves are all almost the same
ax = df[ df['opti'] == 5 ].plot(x='dim', y='loop RbC', label='RbC -O5', color='red', style= '-.',  ax=ax, linewidth=2)

# benchmark with opti
ax = df[ df['opti'] == 1 ].plot(x='dim', y='loop JKI', label='Loop JKI -O1', color='orange', style= '--', ax=ax, linewidth=2)
ax = df[ df['opti'] == 2 ].plot(x='dim', y='loop JKI', label='Loop JKI -O2', color='green',  style= '--', ax=ax, linewidth=2)
ax = df[ df['opti'] == 5 ].plot(x='dim', y='loop JKI', label='Loop JKI -O5', color='dodgerblue', style= '--', 
                                ax=ax, linewidth=2, logy=True)

plt.legend(ncol=2, prop={'size': 21})

# change properties
ax.set(xlabel = 'matrix dimension $N$', ylabel = 'time $[s]$',
       title = "gfortran: optimization flag [$-O*$]", xlim=[500, 10200], ylim=[10**-3,10**4])
ax.grid('on',  linestyle=':', linewidth=1)
plt.yticks(10.0**np.arange(-3,4,1))
plt.tight_layout()

if do_writeplot:
    fname = fname_plots.format("CPUtimes_with_opti")
    print("writing plot", fname)
    plt.savefig( fname, transparent=True)




sys.exit(0)
