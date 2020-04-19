import matplotlib
matplotlib.use('MacOSX')

import matplotlib.pyplot as plt
import numpy as np
import plotlib as plib

def plot_gammas():
    fig, ax = plib.set_plot_style()

    filename = 'test_gamma.txt'
    plotname = 'gamma-xsecs.pdf'
    ax.set_ylabel(r'$E_\gamma d\sigma/dE_\gamma$ [mb]')
    #ax.set_ylabel(r'$E_\nu d\sigma/dE_\nu$ [mb]')

    ax.set_xscale('log')
    ax.set_xlabel('E [GeV]')
    ax.set_xlim([1, 1e5])
    ax.set_xticks([1,10,1e2,1e3,1e4,1e5])
    ax.set_yscale('log')
    ax.set_ylim([1e-1, 1e3])
        
    E, sigma_1e1_K, sigma_1e2_K, sigma_1e3_K, sigma_1e5_K = np.loadtxt(filename, skiprows=0, usecols=(0,1,2,3,4), unpack=True)

    ax.plot(E, E * sigma_1e1_K, color='tab:brown', label='10 GeV')
    ax.plot(E, E * sigma_1e2_K, color='tab:red', label='100 GeV')
    ax.plot(E, E * sigma_1e3_K, color='tab:green', label='1 TeV')
    ax.plot(E, E * sigma_1e5_K, color='tab:blue', label='100 TeV')

    E, sigma_1e1_A, sigma_1e2_A, sigma_1e3_A, sigma_1e5_A = np.loadtxt(filename, skiprows=0, usecols=(0,5,6,7,8), unpack=True)

    ax.plot(E, E * sigma_1e1_A, color='tab:brown', linestyle=':')
    ax.plot(E, E * sigma_1e2_A, color='tab:red', linestyle=':')
    ax.plot(E, E * sigma_1e3_A, color='tab:green', linestyle=':')
    ax.plot(E, E * sigma_1e5_A, color='tab:blue', linestyle=':')

    E, sigma_1e1_M, sigma_1e2_M, sigma_1e3_M, sigma_1e5_M = np.loadtxt(filename, skiprows=0, usecols=(0,9,10,11,12), unpack=True)

    ax.plot(E, E * sigma_1e1_M, color='tab:brown', linestyle='--')
    ax.plot(E, E * sigma_1e2_M, color='tab:red', linestyle='--')
    ax.plot(E, E * sigma_1e3_M, color='tab:green', linestyle='--')
    ax.plot(E, E * sigma_1e5_M, color='tab:blue', linestyle='--')
    
    E, sigma = np.loadtxt('spec_gam', skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(E, sigma, 'k', linestyle=':')
    
    ax.legend(fontsize=14)
    
    plt.savefig(plotname)
    
# MAiN
plot_gammas()
