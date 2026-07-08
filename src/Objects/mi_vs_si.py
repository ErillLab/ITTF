# -*- coding: utf-8 -*-
"""

Class for the analysis of the statistical representation of encoding strategies
for sequence motifs. The code allows

"""

import numpy as np
import itertools
import random
import pandas as pd


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
        Returns the percent SI of all the motifs selecting `n` target sequences.
        This provides the data for the EXACT probability ditribution over %SI.
        
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
        Returns the percent SI of all randomly sampled motifs selecting `n` target sequences.
        This provides the data for the ESTIMATED probability ditribution over %SI.
        
        `sample_size`: size of the sample used to estimate the probability distribution.
                
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

