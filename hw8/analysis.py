# ========================================
#  QIC  >>  assignment 8
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 24 December 2022
# ========================================



# %% import common resources ---------------------------------------

from hw8_common import *

do_writeplot = True



# %% read data from RG & DMRG

rg5 = read_rgobj('data/rg5-30.hw8')
rg6 = read_rgobj('data/rg6-30.hw8')
rg7 = read_rgobj('data/rg7-15.hw8')
rg_data = [ rg5, rg6, rg7 ]


dmrg3b = read_dmrgobj('data/dmrg3-30-big.hw8')
dmrg4  = read_dmrgobj('data/dmrg4-30.hw8')




# %% ======================================================
#  PLOT | plot the spectrum energy of ising model
# =========================================================

plot_spectrum( [ rg6, dmrg3b ],
               labels = ['RG','iDMRG'],
               lw = 3, cl = ['tab:blue', '#FFBB00']
).savefig( fname_plots.format('ground'), transparent=True)  



# %% ======================================================
#  PLOT | iterations required to converge at 1e-11
# =========================================================

plot_iteration_to_convergence( rg_data,
    labels = ['RG $(N=5)$', 'RG $(N=6)$', 'RG $(N=7)$'] ,
    lw = 2
).savefig( fname_plots.format('RG_iter'), transparent=True)


# %% ======================================================
#  PLOT | convergence
# =========================================================

plot_hist( rg5, [ 0, 10, 15, 30 ], cm = 'winter'
).savefig( fname_plots.format('RG_hist'), transparent=True)


# %%

plot_hist( dmrg3b, [ 0, 10, 15, 30 ], cm = 'autumn'
).savefig( fname_plots.format('DMRG_hist'), transparent=True)



# %%
