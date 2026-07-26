# -*- coding: utf-8 -*-
"""

Code to generate Fig. 4C, S7 and S8.

Analysis of the statistical representation of encoding strategies for sequence motifs.


- Fig. 4C:
  Heatmap showing the null probability distribution of %SI, for L=6, annotating
  the steady-state values found in Fig. 4B using special symbols. Only pairwise
  correlations (MI) between adjacent positions are considered.

- Fig. S7:
  Heatmap showing the null probability distribution of %SI, for L values ranging
  from 2 to 8. Only pairwise correlations (MI) between adjacent positions are considered.

- Fig. S8:
  Heatmap showing the null probability distribution of %SI, for L values ranging
  from 2 to 8. All correlations (C) are considered (including non-adjacent and
  non-pairwise, that is, higher-order correlations).

"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from src.Objects.mi_vs_si import MI_VS_SI


def generate_data(L, only_adj_MI=True, use_sampling=True, sample_size=1000):
    ''' Generates the data for the heatmap, that is, the distribution of motifs
    over %SI bins, conditional on total informaion (TI).
    
    `use_sampling`:
        - if True (default), the distribution is estimated via Monte Carlo
          sampling using the specified `sample_size`.
        - if False, the exact distribution is calculated (computationally
          feasible only for small enough L values).
    
    `only_adj_MI`:
        - if True, the percent SI is calculated by considering as
          sources of information only SI and the MI between adjacent motifs:
              percent_SI = 100 * SI / (SI + MI_adj)
        - if False (default), the percent SI is calculated by considering
          all possible correlations:
              percent_SI = 100 * SI / (SI + C)
          where C is the 'total correlation' (Watanabe 1960), also known as
          'multiiformation' (Studený & Vejnarová 1999).
    '''
    # num_Lmers: total number of possibe L-mers
    num_Lmers = 4**L
    # n: number of "selected" (recognized) L-mers
    # n_vals: list of possible values of n that map to integer information values (1, 2, 3, ..., 2L)
    n_vals = [int(num_Lmers/(2**i)) for i in range(1, 2*L + 1)]
    # Loop over `n_vals` to test across information density (from 1 bit to 2L bits)
    df_list = []
    for n in n_vals:
        MVS = MI_VS_SI(L)
        if use_sampling:
            # Estimate distribution via Monte Carlo sampling
            df = MVS.percent_SI_distribution_sample(n, sample_size=sample_size, only_adj_MI=only_adj_MI)
        else:
            # Calculate the exact distribution
            df = MVS.percent_SI_distribution(n, only_adj_MI=only_adj_MI)
        df_list.append(df)
    return pd.concat(df_list, ignore_index=True)

def highlight_cell(ax, x, y, xedges, yedges, color):
    ''' Function to highlight the bin containing a symbol at (x,y). '''
    # Find bin index
    i = np.searchsorted(xedges, x, side='right') - 1
    j = np.searchsorted(yedges, y, side='right') - 1

    # safety check (avoid out-of-range)
    if i < 0 or i >= len(xedges)-1 or j < 0 or j >= len(yedges)-1:
        return
    
    rect = patches.Rectangle((xedges[i], yedges[0]),
                             xedges[i+1] - xedges[i], yedges[-1] - yedges[0],
                             linewidth=0.3, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    rect = patches.Rectangle((xedges[i], yedges[j]),
                             xedges[i+1] - xedges[i], yedges[j+1] - yedges[j],
                             linewidth=1.0, edgecolor=color, facecolor='none')
    ax.add_patch(rect)

def plot_percent_SI_heatmap(final_df, L, ax, annotations=None, vmin=0, vmax=1):
    ''' Makes the heatmap showing the distribution of motifs in terms of %SI. '''
    TI_vals = sorted(set(int(x) for x in final_df['TI']))  # all TI values present in the dataframe
    
    h, xedges, yedges, img = ax.hist2d(final_df['TI'], final_df['percent_SI'],
                                       bins=(len(TI_vals), 8),
                                       range=[[0.5, max(TI_vals) + 0.5], [0, 100]],
                                       vmin=vmin, vmax=vmax)
    # Optional annotations
    if annotations is not None:
        for ann in annotations:
            ax.scatter(ann['x'], ann['y'], marker=ann['marker'], s=ann.get('size', 100),
                       c=ann.get('color', 'black'), edgecolors=ann.get('edgecolors', None),
                       linewidths=ann.get('linewidths', 1), zorder=10)
            if ann.get('highlight', False):
                highlight_cell(ax, ann['x'], ann['y'], xedges, yedges, color='white')
    # Axes
    ax.set_xlabel('TI (bits)')
    ax.set_ylabel('Percent SI')
    ax.set_xticks(TI_vals)
    ax.set_title(f'L = {L}')
    return img

def add_colorbar_in_empty_ax(fig, ax, img, sample_size):
    ''' Colorbar to be shared in multi-pane Fig. S7 and S8. '''
    ax.axis("off")
    norm = mpl.colors.Normalize(vmin=0, vmax=sample_size)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=img.get_cmap())
    sm.set_array([])
    cax = ax.inset_axes([0, 0, 0.1, 1])
    cbar = fig.colorbar(sm, cax=cax)
    ticks = [0, 0.25, 0.5, 0.75, 1.0]
    tick_vals = [t * sample_size for t in ticks]
    cbar.set_ticks(tick_vals)
    cbar.set_ticklabels([f"{t:.2f}" for t in ticks])
    cbar.set_label("p(%SI | TI)", fontsize=14, rotation=270, labelpad=15)

def multipane_sweep_over_L(filename, only_adj_MI, use_sampling, sample_size, annotations=None):
    ''' Creates a multipane plot with one heatmap per value of L (2 to 8 bp). '''
    
    L_values = list(range(2, 9))
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), dpi=600)
    axes = axes.flatten()
    img = None
    
    for ax, L in zip(axes[:-1], L_values):
        # Data
        print(f'Generating data for L={L} ...')
        df = generate_data(L=L, only_adj_MI=only_adj_MI, use_sampling=use_sampling, sample_size=sample_size)
        # Annotations to be passed to 'plot_percent_SI_heatmap'
        ann_to_pass = None
        if L == 3:
            if annotations is not None:
                # Annotate the position of the standard genetic code
                ann_to_pass = annotations
        # Heatmap
        img = plot_percent_SI_heatmap(df, L=L, ax=ax, annotations=ann_to_pass, vmin=0, vmax=sample_size)
    
    # last axis reserved for colorbar
    add_colorbar_in_empty_ax(fig, axes[-1], img, sample_size)
    
    plt.tight_layout()
    plt.savefig(filename, bbox_inches="tight")
    plt.show()

def make_figure_S7():
    multipane_sweep_over_L("Fig_S7.png", only_adj_MI=True, use_sampling=True, sample_size=1000)

def make_figure_S8():
    multipane_sweep_over_L("Fig_S8.png", only_adj_MI=False, use_sampling=True, sample_size=1000,
                           annotations=FIG_S8_L3_ANNOTATIONS)
 
def make_figure_4C():
    ''' Produces Fig. S8. '''
    L = 6
    sample_size = 1000
    
    # Data
    df = generate_data(L=L, only_adj_MI=True, use_sampling=True, sample_size=sample_size)
    
    # Heatmap
    fig, ax = plt.subplots(figsize=(5.5, 4), dpi=600)
    img = plot_percent_SI_heatmap(df, L=L, ax=ax, annotations=FIG_4C_ANNOTATIONS, vmin=0, vmax=sample_size)
    
    # Colorbar
    norm = mpl.colors.Normalize(vmin=0, vmax=sample_size)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=img.get_cmap())
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    ticks = [0, 0.25, 0.5, 0.75, 1.0]
    cbar.set_ticks([t * sample_size for t in ticks])
    cbar.set_ticklabels([f"{t:.2f}" for t in ticks])
    cbar.set_label("p(%SI | TI)", rotation=270, labelpad=20)

    plt.tight_layout()
    plt.savefig("Fig_4C.png", bbox_inches="tight")
    plt.show()


# Figure-specific annotations

star  = {'x': 6.47, 'y': 75.92638228212209, 'marker': '*', 'size': 300,
         'color': 'black', 'linewidths': 0.3, 'highlight': True, 'edgecolors': 'white'}
cross = {'x': 1.64, 'y': 26.51257056268843, 'marker': 'x', 'size': 150,
         'color': 'black', 'linewidths': 2.0, 'highlight': True, 'edgecolors': None}
diam  = {'x': 4.39, 'y': 97.39552627908698, 'marker': 'd', 'size': 150,
         'color': 'red',   'linewidths': 0.3, 'highlight': True, 'edgecolors': 'white'}

FIG_4C_ANNOTATIONS = [star, cross]
FIG_S8_L3_ANNOTATIONS = [diam]


def main():
    make_figure_4C()
    make_figure_S7()
    make_figure_S8()


if __name__ == "__main__":
    main()











