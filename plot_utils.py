import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

np.random.seed(1738)

plt.rcParams.update({
    'font.family':'serif',
    'font.size':14.0,
    'axes.labelsize': 'large',
    'xtick.labelsize':'large',
    'ytick.labelsize':'large',
    'axes.titlesize':'large',
    'legend.fontsize': 'large',
    'lines.linewidth':1.8,
    'patch.linewidth':1.8,})

colors = list(matplotlib.colors.TABLEAU_COLORS.keys())
colors_bin = ['darkslategray', 'forestgreen', 'goldenrod', 'purple']
colors_stell_bin = ['tab:blue', 'tab:red']

cmap_blue = LinearSegmentedColormap.from_list('gray_to_blue', [(0.8, 0.8, 0.8), (31/255, 119/255, 180/255)], N=256)
cmap_orange = LinearSegmentedColormap.from_list('gray_to_orange', [(0.8, 0.8, 0.8), (255/255, 127/255, 14/255)], N=256)
cmap_green = LinearSegmentedColormap.from_list('gray_to_orange', [(0.8, 0.8, 0.8), (44/255, 160/255, 44/255)], N=256)
cmap_red = LinearSegmentedColormap.from_list('gray_to_orange', [(0.8, 0.8, 0.8), (214/255, 39/255, 40/255)], N=256)
colormaps = [cmap_blue, cmap_orange, cmap_green, cmap_red]

colormaps_bin = []
for i, color in enumerate(colors_bin):
    cmap = LinearSegmentedColormap.from_list(f'gray_to_{color}', [(0.8, 0.8, 0.8), color], N=256)
    colormaps_bin.append(cmap)

z_params = {
    'xlim': (0, 3),
    'ylim': (0, 7),
    'xlabel': r'$z$',
    'ylabel': 'n'
}

stell_params = {
    'xlim': (5, 13),
    'ylim': (0, 0.7),
    'xlabel': r'log$_{10}$(M$_*$/M$_{\odot}$)',
    'ylabel': 'n'
}

ssfr_params = {
    'xlim': (-20, -7.5),
    'ylim': (0, 1.6),
    'xlabel': r'log$_{10}$(SSFR / yr$^{-1}$)',
    'ylabel': 'n'
}


def axis_settings(ax, params):
    ax.set_xlim(*params['xlim'])
    ax.set_ylim(*params['ylim'])
    ax.set_xlabel(params['xlabel'])
    ax.set_ylabel(params['ylabel'])
    

def plot_vline(ax, x, ymin=-30, ymax=30, ls='--', lw=1.2, color='black', label=None):
    ax.vlines(x, ymin, ymax, color=color, label=label, ls=ls, lw=lw)

def plot_hline(ax, y, xmin=-30, xmax=30, ls='--', lw=1.2, color='black', label=None):
    ax.hlines(y, xmin, xmax, color=color, label=label, ls=ls, lw=lw)


def flux2mag(flux, zero_pt=30):
    # Convert flux to magnitude
    return zero_pt - 2.5 * np.log10(flux) 
    