import matplotlib
matplotlib.use('MacOSX')

import matplotlib.pyplot as plt
import numpy as np
import plotlib as plib

def plot_gammas():
    fig, ax = plib.set_plot_style()

    plotname = 'gamma-xsecs-test.pdf'
    ax.set_ylabel(r'$E_\gamma d\sigma/dE_\gamma$ [mb]')
    #ax.set_ylabel(r'$E_\nu d\sigma/dE_\nu$ [mb]')

    ax.set_xscale('log')
    ax.set_xlabel('E [GeV]')
    ax.set_xlim([.1, 1e3])
    #ax.set_xticks([1,10,1e2,1e3]),1e4,1e5])
    ax.set_yscale('log')
    ax.set_ylim([1e-1, 1e3])
    
    E, s1e1, s1e2, s1e3, s1e5 = np.loadtxt('specgamma-kamae.txt', skiprows=1, usecols=(0,1,2,3,4), unpack=True)

    #ax.plot(E, E * s1e1, color='tab:brown', label='10 GeV')
    #ax.plot(E, E * s1e2, color='tab:red', label='100 GeV')
    ax.plot(E, E * s1e3, color='tab:green', label='Kamae') # label='1 TeV',
    #ax.plot(E, E * s1e5, color='tab:blue', label='100 TeV')
    
    #E, s1e1, s1e2, s1e3, s1e5 = np.loadtxt('specnu-aafrag.txt', skiprows=1, usecols=(0,1,2,3,4), unpack=True)

    #ax.plot(E, E * s1e1, color='tab:brown', linestyle=':')
    #ax.plot(E, E * s1e2, color='tab:red', linestyle=':')
    #ax.plot(E, E * s1e3, color='tab:green', linestyle=':')
    #ax.plot(E, E * s1e5, color='tab:blue', linestyle=':')
    
    E, s1e1, s1e2, s1e3, s1e5 = np.loadtxt('specgamma-kelner.txt', skiprows=1, usecols=(0,1,2,3,4), unpack=True)

    #ax.plot(E, E * s1e1, color='tab:brown', linestyle='--')
    #ax.plot(E, E * s1e2, color='tab:red', linestyle='--')
    ax.plot(E, E * s1e3, color='tab:green', linestyle='--', label='Kelner \& Aharonian')
    #ax.plot(E, E * s1e5, color='tab:blue', linestyle='--')
    
    E, s1e3 = np.loadtxt('spec_gam', skiprows=0, usecols=(0,1), unpack=True)
    
    ax.plot(E, s1e3, color='tab:green', linestyle=':', label='AAfrag100')

    ax.legend(fontsize=16)
    
    plt.savefig(plotname)
    
def plot_ppgam():
    fig, ax = plib.set_plot_style()

    filename = 'specgamma-ppgam.txt'
    plotname = 'ppgam-xsecs.pdf'
    ax.set_ylabel(r'$E_\gamma d\sigma/dE_\gamma$ [mb]')

    ax.set_xscale('log')
    ax.set_xlabel('E [GeV]')
    ax.set_xlim([1, 1e5])
    ax.set_xticks([1,10,1e2,1e3,1e4,1e5])
    ax.set_yscale('log')
    ax.set_ylim([1e-1, 1e3])

    E, s1e1, s1e2, s1e3, s1e5 = np.loadtxt(filename, skiprows=1, usecols=(0,1,2,3,4), unpack=True)

    ax.plot(E, E * s1e1, color='tab:brown', label='10 GeV')
    ax.plot(E, E * s1e2, color='tab:red', label='100 GeV')
    ax.plot(E, E * s1e3, color='tab:green', label='1 TeV')
    ax.plot(E, E * s1e5, color='tab:blue', label='100 TeV')

    E, s1e1, s1e2, s1e3, s1e5 = np.loadtxt(filename, skiprows=1, usecols=(0,5,6,7,8), unpack=True)

    ax.plot(E, E * s1e1, color='tab:brown', linestyle=':')
    ax.plot(E, E * s1e2, color='tab:red', linestyle=':')
    ax.plot(E, E * s1e3, color='tab:green', linestyle=':')
    ax.plot(E, E * s1e5, color='tab:blue', linestyle=':')

    E, s1e1, s1e2, s1e3, s1e5 = np.loadtxt(filename, skiprows=1, usecols=(0,9,10,11,12), unpack=True)

    ax.plot(E, E * s1e1, color='tab:brown', linestyle='--')
    ax.plot(E, E * s1e2, color='tab:red', linestyle='--')
    ax.plot(E, E * s1e3, color='tab:green', linestyle='--')
    ax.plot(E, E * s1e5, color='tab:blue', linestyle='--')

    E, s1e1, s1e2, s1e3, s1e5 = np.loadtxt(filename, skiprows=1, usecols=(0,13,14,15,16), unpack=True)

    ax.plot(E, E * s1e1, color='tab:brown', linestyle='-.')
    ax.plot(E, E * s1e2, color='tab:red', linestyle='-.')
    ax.plot(E, E * s1e3, color='tab:green', linestyle='-.')
    ax.plot(E, E * s1e5, color='tab:blue', linestyle='-.')

    ax.legend(fontsize=16)

    plt.savefig(plotname)

# MAiN
plot_gammas()
#plot_ppgam()
