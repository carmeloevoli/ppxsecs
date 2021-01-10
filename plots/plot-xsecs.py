import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use("./.mplstyle")
import numpy as np

def plot_model(ax, filename, color, label):
    E, s1e2, s1e3, s1e5 = np.loadtxt(filename, skiprows=1, usecols=(0,2,3,5), unpack=True)
    ax.plot(E / 1e2, E * s1e2, color=color, linestyle='-', label=label)
    ax.plot(E / 1e3, E * s1e3, color=color, linestyle=':')
    ax.plot(E / 1e5, E * s1e5, color=color, linestyle='--')

def plot_ratio(ax, filegamma, filenu, color, label):
    a, rg1, rg2, rg3 = np.loadtxt(filegamma, skiprows=1, usecols=(0,1,2,3), unpack=True)
    a, rn1, rn2, rn3 = np.loadtxt(filenu, skiprows=1, usecols=(0,1,2,3), unpack=True)
    ax.plot(a, rg1 / rn1, color=color, linestyle='-', label=label)
    ax.plot(a, rg2 / rn2, color=color, linestyle=':')
    #ax.plot(a, rg3 / rn3, color=color, linestyle='--')

def plot_model_primary(ax, filename, color, label):
    E, s1e1, s1e3, s1e5 = np.loadtxt(filename, skiprows=1, usecols=(0,1,3,5), unpack=True)
    ax.plot(E / 1e0, 1e1 * s1e1, color=color, linestyle='-', label=label)
    ax.plot(E / 1e0, 1e3 * s1e3, color=color, linestyle=':')
    ax.plot(E / 1e0, 1e5 * s1e5, color=color, linestyle='--')

def plot_gammas(doGamma):
    fig = plt.figure(figsize=(9.5,7.5))
    ax = fig.add_subplot(111)

    if doGamma:
        plotname = 'gamma-xsecs-test.pdf'
        ax.set_ylabel(r'$E_\gamma d\sigma/dE_\gamma$ [mb]')
        ax.set_xlabel(r'$E_\gamma / E_p$')
        ax.set_ylim([0, 200])
    else:
        plotname = 'nu-xsecs-test.pdf'
        ax.set_ylabel(r'$E_\nu d\sigma/dE_\nu$ [mb]')
        ax.set_xlabel(r'$E_\nu / E_p$')
        ax.set_ylim([0, 500])

    ax.set_xscale('log')
    ax.set_xlim([1e-5, 1e0])
    ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    #ax.set_yscale('log')
    
    if doGamma:
        plot_model(ax, '../build/specgamma-kamae.txt', 'tab:red', 'Kamae+06')
        plot_model(ax, '../build/specgamma-aafrag.txt', 'tab:green', 'AAfrag')
        plot_model(ax, '../build/specgamma-kelner.txt', 'tab:orange', 'Kelner')
        #plot_model(ax, '../build/specgamma-ppgam-GEANT4.txt', 'tab:cyan', 'AAfrag')
        plot_model(ax, '../build/specgamma-ppgam-PYTHIA8.txt', 'tab:cyan', 'ppgam-PYTHIA8')
        plot_model(ax, '../build/specgamma-ppgam-QGSJET.txt', 'tab:purple', 'ppgam-QGSJET')
        plot_model(ax, '../build/specgamma-ppgam-SIBYLL.txt', 'tab:blue', 'ppgam-SIBYLL')
    else:
        plot_model(ax, '../build/specnu-kamae.txt', 'tab:red', 'Kamae+06')
        plot_model(ax, '../build/specnu-aafrag.txt', 'tab:green', 'AAfrag')
        plot_model(ax, '../build/specnu-kelner.txt', 'tab:orange', 'Kelner')
        plot_model(ax, '../build/specnu-ppgam-PYTHIA8.txt', 'tab:cyan', 'ppgam-PYTHIA8')
        plot_model(ax, '../build/specnu-ppgam-QGSJET.txt', 'tab:purple', 'ppgam-QGSJET')
        plot_model(ax, '../build/specnu-ppgam-SIBYLL.txt', 'tab:blue', 'ppgam-SIBYLL')

    ax.legend(fontsize=10)
    plt.savefig(plotname, dpi=300)
    
def plot_gammanu_ratio():
    fig = plt.figure(figsize=(9.5,7.5))
    ax = fig.add_subplot(111)

    plotname = 'gammanu-xsecs-test.pdf'
    ax.set_ylabel(r'$\gamma/\nu$')
    ax.set_xlabel(r'$\alpha$')
#    ax.set_ylim([0, 500])
    ax.set_xlim([2,3])
    
    plot_ratio(ax, '../build/rategamma-kamae.txt', '../build/ratenu-kamae.txt', 'tab:red', 'Kamae+06')
    plot_ratio(ax, '../build/rategamma-aafrag.txt', '../build/ratenu-aafrag.txt', 'tab:green', 'AAfrag')
    plot_ratio(ax, '../build/rategamma-kelner.txt', '../build/ratenu-kelner.txt', 'tab:orange', 'Kelner')

    plt.savefig(plotname, dpi=300)

# MAiN
plot_gammas(True)
plot_gammas(False)
plot_gammanu_ratio()
