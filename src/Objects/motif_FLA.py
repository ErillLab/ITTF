# -*- coding: utf-8 -*-
"""



"""


import math
import itertools
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import dipwmsearch as ds

import sys
sys.path.append('../helper_functions')
from helperFunctions import create_dipwm_object



class MotifFLA():
    
    def __init__(self, model):
        self.model = model
        self.length = model.length_bs
        self.kmer_to_score = None
        self.global_opt = None
        self.set_kmer_to_score()
        self.set_global_opt()
        self._bases = ['A','C','G','T']
    
    def get_all_kmers(self):
        '''
        Returns a list of all possible DNA k-mers (as strings), where k is the
        length of the motif model.
        '''
        all_kmers = []
        for kmer in itertools.product('ACGT', repeat=self.length):
            all_kmers.append(''.join(kmer))
        return all_kmers
    
    def set_kmer_to_score(self):
        '''
        Sets the `kmer_to_score` attribute to a dictionary where the keys are
        the k-mers (str) and the values are scores (float).
        '''
        all_kmers = self.get_all_kmers()
        all_scores = [self.model.score(kmer) for kmer in all_kmers]  # XXX
        self.kmer_to_score = dict(zip(all_kmers, all_scores))
    
    def set_global_opt(self):
        '''
        Sets the `global_opt` attribute to a list of k-mers. If the global
        optimum is unique, the list will contain a single k-mer (str).
        Otherwise, it will contain all the k-mers with maximum score.
        '''
        max_score = - math.inf
        global_opt_list = []
        for kmer in self.kmer_to_score.keys():
            score = self.kmer_to_score[kmer]
            if score > max_score:
                max_score = score
                # Reset list of global optima
                global_opt_list = []
                global_opt_list.append(kmer)
            elif score == max_score:
                # Add to list of global optima
                global_opt_list.append(kmer)
        self.global_opt = global_opt_list
    
    def hamming_d(self, seq1, seq2):
        ''' Returns the Hamming distance between the two input strings. '''
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    
    def gdpdc(self, corr='pearson'):
        '''
        Computes the 'Genotype Distance - Phenotype Distance Correlation' (GDPDC).
        For each pair of kmers, it computes:
            (1) the Hamming distance [genotype distance]
            (2) the difference in score (in absolute value) [phenotype distance]
        Then, it returns the correlation (Pearson or Spearman) between (1) and (2).
        '''
        all_kmers = self.kmer_to_score.keys()
        distances = []
        delta_scores = []
        for i in range(len(all_kmers)):
            kmer1 = all_kmers[i]
            for j in range(i+1, len(all_kmers)):
                kmer2 = all_kmers[j]
                # Hamming distance between the two k-mers
                distances.append(self.hamming_d(kmer1, kmer2))
                # Delta-score
                delta_scores.append(abs(self.kmer_to_score[kmer1] - self.kmer_to_score[kmer2]))
        if corr == 'pearson':
            r, pval = pearsonr(distances, delta_scores)
        elif corr == 'spearman':
            r, pval = spearmanr(distances, delta_scores)
        else:
            raise ValueError("corr must be 'pearson' or 'spearman.")
        return r, pval
    
    def fdc(self, corr='pearson'):
        '''
        Computes the 'Fitness - Distance Correlation' (FDC).
        For each kmer, it computes:
            (1) the Hamming distance from the global optimum
            (2) the fitness (the score)
        Then, it returns the correlation (Pearson or Spearman) between (1) and (2).
        '''
        f_list = []
        d_list = []
        for kmer, score in self.kmer_to_score.items():
            # Fitness is the score
            f_list.append(score)
            # Distance from (closest) global optimum
            d_list.append(min([self.hamming_d(kmer, opt)] for opt in self.global_opt))
        if corr == 'pearson':
            r, pval = pearsonr(f_list, d_list)
        elif corr == 'spearman':
            r, pval = spearmanr(f_list, d_list)
        else:
            raise ValueError("corr must be 'pearson' or 'spearman.")
        return r, pval
    
    def no(self):
        '''
        Returns the Number of Optima (NO) in the fitness landscape (both local
        and global), as well as the optima (as a list of sets).
        '''
        optima = []
        all_kmers = self.kmer_to_score.keys()
        for kmer in all_kmers:
            score = self.kmer_to_score[kmer]
            neighbors = [km for km in all_kmers if self.hamming_d(kmer, km)==1]
            # Candidate optimum
            optimum = [kmer]
            for n in neighbors:
                n_score = self.kmer_to_score[n]
                if n_score > score:
                    optimum = False
                    break
                elif n_score == score:
                    optimum.append(n_score)
            if optimum:
                # True optimum. Turn it into a set
                optimum = set(optimum)
                # Add to list if not found yet
                if optimum not in optima:
                    optima.append(optimum)
        return len(optima), optima
    
    def fc(self):
        '''
        Analyzes the fitness landscape with the Fitness Cloud (FC) approach.
        For each kmer it computes:
            (1) Its fitness (score)
            (2) The fitness of all its neighbors (called the 'cloud')
        It compares (1) with (2) [scatter-plot], as well as the average of (2)
        [line-plot]. This plot can be used to understand the evolvability,
        because it shows the likelihood of adaptive variation.
        '''
        all_kmers = self.kmer_to_score.keys()
        f = []
        c = []
        c_avg = []
        for kmer in all_kmers:
            # k-mer fitness
            f.append(self.kmer_to_score[kmer])
            neighbors = [km for km in all_kmers if self.hamming_d(kmer, km)==1]
            # Cloud (fitness of the neighbors)
            cloud = [self.kmer_to_score[n] for n in neighbors]
            c.append(cloud)
            # Cloud average
            c_avg.append(sum(cloud)/len(cloud))
        
        # Plot fitness cloud
        
        # Plot line for cloud average
        plt.plot(f, c_avg, color='red')
        # Plot all the clouds (all the datapoints)
        for i in range(len(f)):
            fitness = f[i]
            cloud = c[i]
            plt.scatter([fitness]*len(cloud), cloud, alpha=0.5)
        plt.show()














