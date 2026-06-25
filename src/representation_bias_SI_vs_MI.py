# -*- coding: utf-8 -*-
"""

Analysis of the statistical representation of encoding strategies for sequence motifs.

"""

import numpy as np
import itertools
import random
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches



class MI_VS_SI():
    '''
    Class to study the representation of encoding strategies with different %SI.
    In other words, how many ways there are to encode information using a given
    encoding strategy (e.g., using only SI, etc.). Statistically, what encoding
    strategy is more "entropic", therefore expected to be frequent when randomly
    sampling strategy space (concept of "representation bias").
    '''
    
    def __init__(self, L):
        '''
        `L`: length of the motif
        '''
        self.L = L
        self.num_Lmers = 4 ** L
        self.all_Lmers = [''.join(seq) for seq in itertools.product('TCAG', repeat=self.L)]
    
    def percent_SI_distribution(self, n, only_adj_MI=False):
        '''
        Returns the exact percent SI.
        
        `n`: number of chosen sequences (recognized by the TF)
        
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
        # TI: total information content (in bits)
        self.TI = np.log2(self.num_Lmers/n)
        # Iterate through all possible choices and compute percent SI
        percent_SI_list = []
        for seqs in itertools.combinations(self.all_Lmers, n):
            # Information encoded through SI
            SI = self._compute_SI(seqs)
            # Define the total as TI or just the sum of SI and adj_MI
            if only_adj_MI:
                # Information encoded through MI between adjacent positions
                adj_MI = self._compute_adj_MI(seqs)
                # total info considered
                total = SI + adj_MI
            else:
                total = self.TI
            # Store the %SI
            percent_SI_list.append(100 * SI / total)
        return pd.DataFrame({'L': self.L, 'n': n, 'TI': self.TI, 'percent_SI': percent_SI_list})
    
    def percent_SI_distribution_sample(self, n, sample_size=10000, only_adj_MI=False):
        '''
        Returns the estimated percent SI, based on a sample of specified size.
        
        `n`: number of chosen sequences (recognized by the TF)
        
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
        # TI: total information content (in bits)
        self.TI = np.log2(self.num_Lmers/n)
        # Iterate through all possible choices and compute percent SI
        percent_SI_list = []
        for i in range(sample_size):
            seqs = random.sample(self.all_Lmers, n)
            # Information encoded through SI
            SI = self._compute_SI(seqs)
            # Define the total as TI or just the sum of SI and adj_MI
            if only_adj_MI:
                # Information encoded through MI between adjacent positions
                adj_MI = self._compute_adj_MI(seqs)
                # total info considered
                total = SI + adj_MI
            else:
                total = self.TI
            # Store the %SI
            percent_SI_list.append(100 * SI / total)
        return pd.DataFrame({'L': self.L, 'n': n, 'TI': self.TI, 'percent_SI': percent_SI_list})
    
    def _compute_SI(self, seqs):
        # Variable to store SI
        SI = 0
        n = len(seqs)
        for col in range(len(seqs[0])):
            # Shannon entropy
            H = 0
            # Iterate through rows to get all bases at a column
            G_count = 0
            C_count = 0
            A_count = 0
            T_count = 0
            for row in range(n):
                nucl = seqs[row][col]
                if nucl == 'G':
                    G_count += 1
                elif nucl == 'C':
                    C_count += 1
                elif nucl == 'A':
                    A_count += 1
                elif nucl == 'T':
                    T_count += 1
                else:
                    raise ValueError("Unknown char")
            # Calculate frequencies
            G_count /= n
            C_count /= n
            A_count /= n
            T_count /= n
            # Calculate uncertainty of the column
            if G_count != 0:
                H -= G_count * np.log2(G_count)
            if C_count != 0:
                H -= C_count * np.log2(C_count)
            if A_count != 0:
                H -= A_count * np.log2(A_count)
            if T_count != 0:
                H -= T_count * np.log2(T_count)
            # Add Positional Information (SI) of the column to the cumulative SI
            SI += 2 - H
        return SI
    
    def _compute_adj_MI(self, seqs):
        adj_MI = 0
        # iterate over all the L-1 adjacent pairs
        for i in range(self.L - 1):
            dimers = [seq[i:i+2] for seq in seqs]
            dimer_counts = [dimers.count(i) for i in set(dimers)]
            dimer_freqs = [c/len(dimers) for c in dimer_counts]
            # Entropy of the dimers at pos (i,i+1)
            H = 0
            for f in dimer_freqs:
                H -= f * np.log2(f)
            # Total information of the dimers at pos (i,i+1)
            TI_dimers = 4 - H  # Because Hbefore is log2(4^L)=log2(16)=4 bits
            # SI of the dimers at pos (i,i+1)
            SI_dimers = self._compute_SI(dimers)
            # MI between pos i and i+1 (MI of this dimer)
            MI_dimers = TI_dimers - SI_dimers
            # Add MI from this dimer to the overall MI between adjacent positions
            adj_MI += MI_dimers
        return adj_MI


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

def make_figure_S7(use_sampling=True):
    ''' Produces Fig. S7. '''
    L_values = list(range(2, 9))
    sample_size = 100
    
    fig, axes = plt.subplots(2, 4, figsize=(12, 6), dpi=600)
    axes = axes.flatten()
    img = None
    
    for ax, L in zip(axes[:-1], L_values):
        # Data
        print(f'Generating data for L={L} ...')
        df = generate_data(L=L, only_adj_MI=True, use_sampling=use_sampling, sample_size=sample_size)
        # Heatmap
        img = plot_percent_SI_heatmap(df, L=L, ax=ax, vmin=0, vmax=sample_size)
    
    # last axis reserved for colorbar
    add_colorbar_in_empty_ax(fig, axes[-1], img, sample_size)

    plt.tight_layout()
    plt.savefig("Figure_S7.png", bbox_inches="tight")
    plt.show()

def make_figure_S8(use_sampling=True):
    ''' Produces Fig. S8. '''
    L_values = list(range(2, 9))
    sample_size = 100
    
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), dpi=600)
    axes = axes.flatten()
    img = None
    
    for ax, L in zip(axes[:-1], L_values):
        # Data
        print(f'Generating data for L={L} ...')
        df = generate_data(L=L, only_adj_MI=False, use_sampling=use_sampling, sample_size=sample_size)
        # Annotate the position of the standard genetic code
        annotations = None
        if L == 3:
            annotations = FIG_S8_L3_ANNOTATIONS
        # Heatmap
        img = plot_percent_SI_heatmap(df, L=L, ax=ax, annotations=annotations, vmin=0, vmax=sample_size)
    
    # last axis reserved for colorbar
    add_colorbar_in_empty_ax(fig, axes[-1], img, sample_size)
    
    plt.tight_layout()
    plt.savefig("Figure_S8.png", bbox_inches="tight")
    plt.show()

def make_figure_4C(use_sampling=True):
    ''' Produces Fig. S8. '''
    L = 6
    sample_size = 1000
    
    # Data
    df = generate_data(L=L, only_adj_MI=True, use_sampling=use_sampling, sample_size=sample_size)
    
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


# Generate Figures
# ----------------

# Figure-specific annotations

star  = {'x': 6.47, 'y': 75.92638228212209, 'marker': '*', 'size': 300,
         'color': 'black', 'linewidths': 0.3, 'highlight': True, 'edgecolors': 'white'}
cross = {'x': 1.64, 'y': 26.51257056268843, 'marker': 'x', 'size': 150,
         'color': 'black', 'linewidths': 2.0, 'highlight': True, 'edgecolors': None}
diam  = {'x': 4.25, 'y': 97, 'marker': 'd', 'size': 150,
         'color': 'red',   'linewidths': 0.3, 'highlight': True, 'edgecolors': 'white'}

FIG_4C_ANNOTATIONS = [star, cross]
FIG_S8_L3_ANNOTATIONS = [diam]


make_figure_4C()
make_figure_S7()
make_figure_S8()











