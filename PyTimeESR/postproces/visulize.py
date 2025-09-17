import sys

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

from .rabi import expand_rates

def plot_rates(fn: str, n: int = 0, globalmax: bool = False,
               levels: list = ['g','e', '0', '2'],
               remove_conserve: bool = True,
               magnitue: bool = False):
    """Plots all terms of Gamma_n. 

    Parameters:
    -----------
    fn (str):   
        File from which to load data 
    n (int, optional):
        Fourier component, default is 0.
    globalmax (bool, optional):
        if True all subplots share the same colormap, default False.
    levels (list, optional):
        list of levels sorted from lowest to highes energy, for tick labels of the plots.
    remove_conserve (bool, optional):
        remove indexes that conserve electron number on one side of Gamma, which must be zero in Sequential theory
    magnitude (bool, optional):
        if True, will only plot magnitudes not sign of gamma components
    
    Returns:
    --------
    fig (plt.figure):
        Figure with the rates ploted
    ax (list): 
        list of axes.
    """
    N = len(levels)
    cmap = 'Blues' if magnitue else  'seismic'

    fig, ax = plt.subplots(2,2, figsize=(12,12))
    fig.tight_layout()
    ax= ax.flatten()

    titles = [r'$\Re{\Gamma_{uv,lj,L}}$',
          r'$\Im{\Gamma_{uv,lj,L}}$',
          r'$\Re{\Gamma_{uv,lj,R}}$',
          r'$\Im{\Gamma_{uv,lj,R}}$',
        ]

    trel_cords = [[[0,2],[2,0]], [[1,2],[2,1]], 
              [[0,3],[3,0]], [[1,3],[3,1]]]
    rabi_cords = [[[0,2],[2,1]], [[1,2],[2,0]],
                  [[0,3],[3,1]], [[1,3],[3,0]]]
    curr_cords = [[[l,j],[j,u]] for l in range(4) for j in range(4) for u in range(4)]
    # listing levels that should be zero
    remove_cords = [[0,0], [1,0], [0,1], [1,1], [2,2], [3,2], [2,3], [3,3]]

    # convert into N**2 by N**2 shape
    trel_cords = np.array([[c[0]+N*c[1] for c in cords] for cords in trel_cords]).T
    rabi_cords = np.array([[c[0]+N*c[1] for c in cords] for cords in rabi_cords]).T
    curr_cords = np.array([[c[0]+N*c[1] for c in cords] for cords in curr_cords]).T
    remove_cords = [c[0]+N*c[1] for c in remove_cords] # converting to 16 by 16
    
    if not remove_conserve:
        # no indexes are removed
        remove_cords = []
    
    # remove indexes
    use_cords = [i for i in range(16) if i not in remove_cords]
    
    # labels for x and y ticks
    labels = np.array([levels[i]+levels[j] for i in range(4) for j in range(4)])
    labels = labels[use_cords]

    # Since we generate the current coordinates automatically some make no sense
    curr_cords = np.array([c for c in curr_cords.T if not c[0] in remove_cords and not c[1] in remove_cords]).T

    # convert coordinates to acount for removin some of them
    rabi_cords = np.array([[c - np.sum(remove_cords<c) for c in rabi_cords[i]] for i in range(2)])
    trel_cords = np.array([[c - np.sum(remove_cords<c) for c in trel_cords[i]] for i in range(2)])
    curr_cords = np.array([[c - np.sum(remove_cords<c) for c in curr_cords[i]] for i in range(2)])
    
    axsize = len(labels) # number of points on the grid
    
    # marker size for highliting tau and omega contributions
    markersize = (3/axsize*fig.dpi)**2 
    
    G = np.loadtxt(fn, skiprows=1)
    Gn = G[G[:,4]==n]
    
    maxmag = np.max(np.abs(Gn[:,5:]))
    minmag = 0. if magnitue else -maxmag

    Gn = np.abs(Gn) if magnitue else Gn
    
    for i in range(4):
        # convert G into a 16x16 coordinace (8x8 if reduced)
        Gni = expand_rates(Gn, N)[:,:,:,:, i]
        Gni16 = Gni.reshape((N**2,N**2))
        Gni16 = Gni16[use_cords][:,use_cords]

        
        # define your own nomralization if asked for
        maxmag = maxmag if globalmax else np.max(abs(Gni))
        minmag = 0. if magnitue else -maxmag
        norm = matplotlib.colors.Normalize(vmin=minmag, vmax=maxmag)
    
        ax[i].set_xticks(np.arange(axsize))
        ax[i].set_yticks(np.arange(axsize))    
        ax[i].set_xticklabels(labels)
        ax[i].set_yticklabels(labels)
    
        ax[i].set_title(titles[i])
        ax[i].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
        ax[i].grid()
        
        ax[i].imshow(Gni16, cmap, norm)

        # Higligth current contributions
        ax[i].scatter(curr_cords[0], curr_cords[1], markersize, marker='s', facecolor='none', 
                  edgecolor='k', linewidth=2.5, label=r'$I$', alpha=.7) 
        
        # Highlight tau contributions
        if i == 0 or i == 2: # only for Re Gamma
            ax[i].scatter(trel_cords[0], trel_cords[1],markersize, marker='s', facecolor='none', 
                  edgecolor='g', linewidth=1.5, label=r'$\tau_{rel}^{-1}$', alpha=.7)

        # Highlight rabi contributions 
        ax[i].scatter(rabi_cords[0], rabi_cords[1], markersize, marker='s', facecolor='none', 
                  edgecolor='orange', linewidth=1.5, label=r'$\Omega$', alpha=.7) 
        
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax[i])

    return fig, ax