# -*- coding: utf-8 -*-
"""

The GeneticCode class constructs a on object that stores the properties of the
standard genetic code or (if randomized) a random code that respects the
structure of the standard genetic code (in terms of number of codons per message),
but with permuted codon-message assignments.

"""


import itertools
import pandas as pd
import math
import numpy as np
import random
import copy



class GeneticCode():
    
    def __init__(self):
        
        self._bases = ['A', 'C', 'G', 'T']
        self.codons = [''.join(triplet) for triplet in itertools.product('TCAG', repeat=3)]
        self._translations = (['Phe'] * 2 + ['Leu'] * 2 + ['Ser'] * 4 + ['Tyr'] * 2 +
                              ['Stop']* 2 + ['Cys'] * 2 + ['Stop']* 1 + ['Trp'] * 1 +
                              ['Leu'] * 4 + ['Pro'] * 4 + ['His'] * 2 + ['Gln'] * 2 +
                              ['Arg'] * 4 + ['Ile'] * 3 + ['Met'] * 1 + ['Thr'] * 4 +
                              ['Asn'] * 2 + ['Lys'] * 2 + ['Ser'] * 2 + ['Arg'] * 2 +
                              ['Val'] * 4 + ['Ala'] * 4 + ['Asp'] * 2 + ['Glu'] * 2 +
                              ['Gly'] * 4)
        self.table = pd.DataFrame({'codons': self.codons, 'translations': self._translations})
        
        self.messages = list(set(self.table['translations']))
        
        # Number of codons per message
        self.num_codons_dict = dict({})
        for message in self.messages:
            self.num_codons_dict[message] = (self.table.translations == message).sum()
        
        self.num_codons_list = list(self.num_codons_dict.values())
        
        # Frequency of each message (amino acids and "stop")
        # Initialized as all equal
        self.message_frequency_dict = dict({})
        for message in self.messages:
            self.message_frequency_dict[message] = 1/len(self.messages)
    
    def set_message_frequencies(self, freq_dict):
        if abs(sum(freq_dict.values()) - 1) > 0.0001:
            raise ValueError("The sum of the amino acid frequencies should be 1")
        if not set(freq_dict.keys()) == set(self.message_frequency_dict.keys()):
            correct_keys = list(self.message_frequency_dict.keys())
            raise ValueError(f"Dictionary keys do not match. They should be {correct_keys}.")
        # Set message frequencies
        for message in freq_dict.keys():
            self.message_frequency_dict[message] = freq_dict[message]
    
    def si_percent(self, exclude_stop=False):
        ''' Returns the overall percentage of information encoded in the form
        of SI (Single-locus Information) in the (genetic) code. '''
        # Messages to consider
        if exclude_stop:
            messages = [m for m in self.messages if m != 'Stop']
        else:
            messages = self.messages
        # Calculate overall proportion of Single-locus Information (SI)
        # print('message , SI , total correlation')
        si_overall = 0
        for message in messages:
            # Codons encoding for `message`
            codons = list(self.table[self.table.translations == message]['codons'])
            # si and mi used to encode `message`
            si, c = self._compute_info(codons)
            # print(message, ',', si, ',', c)
            si_proportion = si / (si + c)
            # Weight by message frequency
            si_overall += (si_proportion * self.message_frequency_dict[message])
        return 100 * si_overall  # return as percentage
    
    def robustness(self, exclude_stop=False):
        ''' Returns the overall mutational robustness of the (genetic) code. '''
        # Messages to consider
        if exclude_stop:
            messages = [m for m in self.messages if m != 'Stop']
        else:
            messages = self.messages
        # Calculate overall mutational robustness
        robustness = 0
        for message in messages:
            # Codons encoding for `message`
            codons = list(self.table[self.table.translations == message]['codons'])
            # Number of correct translations after mutation (initialize at max)
            correct = 9 * len(codons)  # 9 possible mutations per codon
            for codon in codons:
                correct -= self._count_mut_errors(codon)
            correct_proportion = correct / (9 * len(codons))
            # Weight by message frequency
            robustness += (correct_proportion * self.message_frequency_dict[message])
        return 100 * robustness
    
    def randomize_code(self):
        codons = np.copy(self.codons)
        np.random.shuffle(codons)
        self.codons = list(codons)
        self.table['codons'] = self.codons
    
    def swap_two_codons(self):
        codons = self.codons[:]
        i, j = random.sample(range(0, 64), 2)
        codons[i], codons[j] = codons[j], codons[i]
        self.codons = codons
        self.table['codons'] = self.codons
    
    def clone(self):
        return copy.deepcopy(self)
    
    def count_codes():
        '''
        Calculates the number of codes that respect the structure of the standard
        genetic code (in terms of number of codons per message), but with
        permuted codon-message assignments. Therefore, the returned number is the
        size of the set from which the `randomize_code` function is sampling.
        '''
        # Number of possible genetic codes that respect the structure is N = 64! / den
        counts = [2, 3, 2, 2, 4, 2, 2, 2, 4, 1, 3, 2, 2, 4, 2, 6, 4, 4, 6, 6, 1]
        den = 1
        for count in counts:
            den *= math.factorial(count)
        # Logarithm (base-10) of den
        log_den = np.log10(den / 100) + 2
        # Logarithm (base-10) of 64! (64 factorial)
        log_num = 0
        for i in range(1, 65):
            log_num += np.log10(i)
        # Logarithm (base-10) of N
        log_N = log_num - log_den
        # N is 10^log_N
        return 10**log_N
    
    def _count_mut_errors(self, codon):
        # Message encoded by `codon`
        message = self._decode_message(codon)
        # Mutated codons
        mut_codons = self._get_mut_codons(codon)
        # Translate mutated codons
        mut_messages = [self._decode_message(mc) for mc in mut_codons]
        # Number of errors
        return len(set(mut_messages).difference(set([message])))
    
    def _get_mut_codons(self, codon):
        mut_codons = []
        for i in range(3):
            mut_bases = list(set(self._bases).difference(set(codon[i])))
            mut_codons += [codon[:i] + b + codon[i+1:] for b in mut_bases]
        return mut_codons
    
    def _decode_message(self, codon):
        ''' Returns the message that the input `codon` encodes for. '''
        return self.table.loc[self.table['codons'] == codon,'translations'].values[0]
    
    def _compute_SI(self, codons):
        
        num_of_codons = len(codons)
        
        #Variable to store cumulative IC
        SI = 0
        
        for col in range(3):
            # Shannon entropy
            H = 0
            
            # Iterate through rows to get all bases at a column
            G_count = 0
            C_count = 0
            A_count = 0
            T_count = 0
            for row in range(num_of_codons):
                nucl = codons[row][col]
                if nucl == 'G':
                    G_count += 1
                elif nucl == 'C':
                    C_count += 1
                elif nucl == 'A':
                    A_count += 1
                else:
                    T_count += 1
            
            # Calculate frequencies
            G_count /= num_of_codons
            C_count /= num_of_codons
            A_count /= num_of_codons
            T_count /= num_of_codons
            
            # Calculate uncertainty of the column
            if G_count != 0:
                H -= G_count * np.log2(G_count)
            if C_count != 0:
                H -= C_count * np.log2(C_count)
            if A_count != 0:
                H -= A_count * np.log2(A_count)
            if T_count != 0:
                H -= T_count * np.log2(T_count)
            
            # Add Single-locus Information (SI) of the column to the cumulative SI
            SI += 2 - H
        
        return SI
    
    def _compute_TI(self, codons):
        if len(codons) != len(set(codons)):
            raise ValueError('This function assumes each element in codons is unique.')
        return np.log2(64/len(codons))
    
    def _compute_MI(self, codons):
        return self._compute_TI(codons) - self._compute_SI(codons)
    
    def _compute_info(self, codons):
        TI = self._compute_TI(codons)
        SI = self._compute_SI(codons)
        MI = TI - SI
        return (SI, MI)











