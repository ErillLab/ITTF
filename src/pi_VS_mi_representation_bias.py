# -*- coding: utf-8 -*-
"""

Analysis of the statistical representation of encoding strategies for sequence motifs.

"""

import numpy as np
import itertools
import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches



class MI_VS_PI():
    '''
    Class to study the representation of encoding strategies with different %PI.
    In other words, how many ways there are to encode information using a given
    encoding strategy (e.g., using only PI, etc.). Statistically, what encoding
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
    
    def percent_PI_distribution(self, n, only_adj_MI=False):
        '''
        `n`: number of chosen sequences (recognized by the TF)
        '''
        # TI: total information content (in bits)
        self.TI = np.log2(self.num_Lmers/n)
        # Iterate through all possible choices and compute percent PI
        percent_PI_list = []
        for seqs in itertools.combinations(self.all_Lmers, n):
            # Information encoded through PI
            PI = self._compute_PI(seqs)
            # Define the total as TI or just the sum of PI and adj_MI
            if only_adj_MI:
                # Information encoded through MI between adjacent positions
                adj_MI = self._compute_adj_MI(seqs)
                # total info considered
                total = PI + adj_MI
            else:
                total = self.TI
            # Store the %PI
            percent_PI_list.append(100 * PI / total)
            '''
            Insert here some code for neighborhood analysis?
            '''
        return pd.DataFrame({'L': self.L, 'n': n, 'TI': self.TI, 'percent_PI': percent_PI_list})
    
    def percent_PI_distribution_sample(self, n, sample_size=10000, only_adj_MI=False):
        '''
        `n`: number of chosen sequences (recognized by the TF)
        '''
        # TI: total information content (in bits)
        self.TI = np.log2(self.num_Lmers/n)
        # Iterate through all possible choices and compute percent PI
        percent_PI_list = []
        for i in range(sample_size):
            seqs = random.sample(self.all_Lmers, n)
            # Information encoded through PI
            PI = self._compute_PI(seqs)
            # Define the total as TI or just the sum of PI and adj_MI
            if only_adj_MI:
                # Information encoded through MI between adjacent positions
                adj_MI = self._compute_adj_MI(seqs)
                # total info considered
                total = PI + adj_MI
            else:
                total = self.TI
            # Store the %PI
            percent_PI_list.append(100 * PI / total)
            '''
            Insert here some code for neighborhood analysis?
            '''
        return pd.DataFrame({'L': self.L, 'n': n, 'TI': self.TI, 'percent_PI': percent_PI_list})
    
    def _compute_PI(self, seqs):
        # Variable to store PI
        PI = 0
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
            # Add Positional Information (PI) of the column to the cumulative PI
            PI += 2 - H
        return PI
    
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
            # PI of the dimers at pos (i,i+1)
            PI_dimers = self._compute_PI(dimers)
            # MI between pos i and i+1 (MI of this dimer)
            MI_dimers = TI_dimers - PI_dimers
            # Add MI from this dimer to the overall MI between adjacent positions
            adj_MI += MI_dimers
        return adj_MI




# ============================================
# Set value of L and run a param sweep over TI
# ============================================

# Settings
# --------

L = 6  # Motif length
sample_size = 1000  # Sample size for probability estimation
only_adj_MI = True  # Whether to consider only the MI between adjacent position or all possible pairs


# Representation analysis
# -----------------------

# num_Lmers: total number of possibe L-mers
num_Lmers = 4**L
# n: number of "selected" (recognized) L-mers
# n_vals: list of possible values of n that map to integer information values (1, 2, 3, ..., 2L)
n_vals = [int(num_Lmers/(2**i)) for i in range(1, 2*L + 1)]
# Loop over `n_vals` to test across information density (from 1 bit to 2L bits)
df_list = []
for n in n_vals:
    MVP = MI_VS_PI(L)
    df = MVP.percent_PI_distribution_sample(n, sample_size=sample_size, only_adj_MI=only_adj_MI)
    df_list.append(df)
final_df = pd.concat(df_list)

# Plot results
TI_vals = list(set([int(x) for x in final_df['TI']]))  # all TI values present in the dataframe
plt.title('L = ' + str(L))
plt.hist2d(final_df['TI'], final_df['percent_PI'], bins=(len(TI_vals), 8), range=([[0.5,max(TI_vals)+0.5],[0,100]]))
plt.xlabel('TI (bits)')
plt.ylabel('Percent PI')
plt.xticks(TI_vals)
#plt.colorbar()
plt.show()



# ===============================
# Version with symbols for Fig. 4
# ===============================

# Prepare data
TI_vals = sorted(list(set(int(x) for x in final_df['TI'])))

# Symbol coordinates
star = (6.47, 75.92638228212209)
cross = (1.64, 26.512570562688435)

# Plot heatmap
plt.figure(figsize=(5.5,4), dpi=300)
h, xedges, yedges, img = plt.hist2d(final_df['TI'], final_df['percent_PI'],
                                    bins=(len(TI_vals), 8), range=[[0.5, max(TI_vals)+0.5], [0, 100]])
ax = plt.gca()

# Overlay symbols
plt.scatter(star[0],  star[1],  marker='*', s=300, c='black', edgecolors='white', linewidths=0.3, zorder=3)
plt.scatter(cross[0], cross[1], marker='x', s=150, c='black', linewidths=2, zorder=3)


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

# Highlight cells containing the symbols
highlight_cell(ax, star[0], star[1], xedges, yedges, color='white')
highlight_cell(ax, cross[0], cross[1], xedges, yedges, color='white')

# Labels and colorbar
plt.xlabel('TI (bits)')
plt.ylabel('Percent PI')
plt.xticks(TI_vals)
cbar = plt.colorbar()
#cbar.set_label(r'$p(\%PI \mid TI)$', rotation=270, labelpad=20)
cbar.set_label('p($Percent$ $PI$ | $TI$)', rotation=270, labelpad=20)

plt.title(f'L = {L}')
plt.savefig('Fig_4C_v5.png', bbox_inches='tight')
plt.show()


















