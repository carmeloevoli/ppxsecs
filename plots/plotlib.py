import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def set_plot_style():
    #plt.style.use('bmh')
    #print(matplotlib.rcParams.keys())
    matplotlib.rcParams.update({
                               #'axes.grid': True,
                               #'axes.titlesize': 'medium',
                               'font.family': 'serif',
                               'font.serif': 'Palatino', #'Helvetica Neue',
                               'font.size': 33,
                               #'grid.color': 'w',
                               #'grid.linestyle': '-',
                               #'grid.alpha': 0.5,
                               #'grid.linewidth': 1,
                               'legend.frameon': False,
                               'legend.fancybox': False,
                               'legend.fontsize': 20,
                               'legend.numpoints': 1,
                               'legend.loc': 'best',
                               #'legend.framealpha': 0.7,
                               #'legend.handletextpad': 0.1,
                               #'legend.labelspacing': 0.2,
                               'lines.linewidth': 3,
                               'savefig.bbox': 'tight',
                               #'savefig.pad_inches': 0.02,
                               'text.usetex': True,
                               #'text.latex.preamble': r'\usepackage{txfonts}',
                               'xtick.labelsize': 30,
                               'ytick.labelsize': 30,
                               'xtick.direction': 'in',
                               'ytick.direction': 'in',
                               'axes.labelpad': 10,
                               'figure.autolayout': True,
                               })
    fig = plt.figure(figsize=(9.15, 8.7))
    ax = fig.add_subplot(1, 1, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.75)
    ax.minorticks_on()
    ax.tick_params('both', length=15, width=1.5, which='major', pad=10, bottom=True, top=True, left=True, right=True)
    ax.tick_params('both', length=0, width=1.3, which='minor', pad=10, bottom=True, top=True, left=True, right=True)
    return fig, ax

matplotlib.rc("savefig", dpi=200)
#matplotlib.rc("figure", figsize=(9.0, 8.67))

