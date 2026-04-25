import random
import numpy as np
from Bio import SeqIO, motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.motifs import Motif
from math import log2
import json
import itertools

import sys
sys.path.append('../helper_functions')
from src.helperFunctions import create_dipwm_object



'''
    Information in Transcription Factor Binding Sites (TFBS) is typically conceived
    as being conveyed column-wise, measured through Information Content (IC).
    Information, however, can also be conveyed between columns (correlation). Mutual 
    Information measures the information shared between two variables or columns. 
    
    The IT_TF project seeks to explore questions related to evolutionary preferences
    between these two ways of storing information and what this reveals about MI.
    This algorithm attempts to provide insight into this question by 
    modeling evolution in TFBS, fitting populations of organisms 
    based on Total Information (The sum of IC and MI).
Finish documentation for script
???

'''

class MotifModel():
    """ 
    Organism object for IT_TF Genetic Algorithm
    
    The IT_TF Genetic Algorithm seeks to discover properties of MI and IC 
    motifs. This organism object initializes a motif with methods to determine
    total information, mutate the object, and export the organism. Organisms 
    are initialized a Population Object to run the genetic algorithm in order
    to select for a Population of motifs with a certain information encoding strategy
    
    This class possess the following methods:   
    - Constructor to initialize motif organism with passed in attributes
    - Several methods for computing the total information, mutual information,
    and information content of the motif's Transcription Factor Binding Sites
    - Mutate method to alter letters in the binding sites. This is used
    to select for motifs with a certain encoding strategy in the Population Class
    """
    
    def __init__(self, N, L, GC_cont, mut_rate, ID, Pseudocount=None, MI_mode=1, read_in=None):
        """
        Initializes motif organism with arguments initialized as attributes.
        Motif organism features a list of TFBS, which is randomly initialized
        based on bp weighting from GC_cont
        
        Parameters:
            N (int) - Number of binding sites
            L (int) - length of binding sites
            GC_cont (float) - GC content of the organism's genome. 
            Used to weight random initialized of binding sites (0-1)
            Mut_rate (float) - Probability of a base pair mutating (0-1)
            MI_mode (int) - Attribute relavent for computing Mutual Information of TFBS.
            Distance between beginning of pairs for computing MI (1 = overlapping pairs, 2 = independent)
            readin (list) - Optional if reading from a list of TFBS. Should be a list or none 
        
        """
        self.number_bs = N
        self.length_bs = L
        self.GC = GC_cont
        self.mut_rate = mut_rate
        self.ID = ID
        # XXX Pseudocount : unused parameter!
        self.MI_mode = MI_mode
        self.read_in = read_in
        self._alphabet = ["G","C","A","T"]
        
        # Initialize attributes
        self.TFBS = []
        self.IC = 0
        self.MI = 0
        self.TI = 0
        self.cum_joint_entropy = 0
        
        # Set attributes
        self.set_TFBS()
        self.set_IC()
        self.set_MI()
        self.TI = self.IC + self.MI
    
    def set_TFBS(self):
        '''
        Generates random TFBS or reads them from a list or a fasta file,
        depending on `read_in`.
        '''
        
        if self.read_in == None:
            # Generate random TFBS, based on GC content
            tfbs = []
            base_prob = [self.GC/2, self.GC/2, (1-self.GC)/2, (1-self.GC)/2]
            for i in range(self.number_bs):
                tfbs.append(''.join(random.choices(self._alphabet, base_prob, k=self.length_bs)))
            self.TFBS = tfbs
        
        elif type(self.read_in) == list :
            # Read TFBS from list
            self.TFBS = self.read_in
            
        elif type(self.read_in) == str:            
            # Read TFBS from fasta file
            tfbs = []
            fas_iterator = SeqIO.parse(self.read_in, "fasta")
            for bs in fas_iterator:
                tfbs.append(str(bs.seq))
            self.TFBS = tfbs
            
            #Reinitializes variables
            self.number_bs = len(self.TFBS)
            self.length_bs = len(self.TFBS[0])
        
        else:
            raise ValueError("read_in paramter: incorrect type")
    
    def set_IC(self):
        '''
        Organism method used to compute Information Content of the organism's TFBS.
        By calculating Information Content and Mutual Information over many generations
        of the population, the preferences of Transcription Factor Binding Sites towards
        one manner of storing information can be explored
        
        This method calculates column-wise information by iterating through the 
        list of TFBS by column, computing the frequency of each base at a column,
        and calculating shannon entropy given the computed frequencies at each column. 
        
        IC_x = Hg - Hx = 2 - Hx
        '''
        #Variable to store cumulative IC
        self.IC = 0        

        #Iterates through list of TFBS to get column wise frequencies
        
        #Iterates through column to get all bases at each column
        for col in range(self.length_bs):
            
            #Variable to store shannon entropy
            Hx = 0 
            
            #Variables to store counts
            G_count = 0
            C_count = 0
            A_count = 0
            T_count = 0

            #Iterates through row to get all bases at a column
            for row in range(self.number_bs):
                
                #Accesses 2D array elements in a certain column
                bp_TFBS = self.TFBS[row][col]
                
                #Keeps count of the base 
                if bp_TFBS == 'G':
                    G_count += 1
                    
                elif bp_TFBS == 'C':
                    C_count += 1
                
                elif bp_TFBS == 'A':
                    A_count += 1
                    
                else:
                    T_count += 1
                    
            #After counting all the bases in a certain column, calculates frequency
            G_count /= self.number_bs
            C_count /= self.number_bs
            A_count /= self.number_bs
            T_count /= self.number_bs
            
            #Calculates uncertainty of the column after knowing the frequencies of each base
            if G_count != 0:
                Hx -= G_count * np.log2(G_count)
            
            if C_count != 0:
                Hx -= C_count * np.log2(C_count)
            
            if A_count != 0:
                Hx -= A_count * np.log2(A_count)
                
            if T_count != 0:
                Hx -= T_count * np.log2(T_count)
            
            #Calculates Information Content of the column and adds it to cumulative IC variable
            self.IC += 2 - Hx
    
    def set_MI(self):
        '''
        Organism method used to compute Mutual Information of the organism's TFBS.
        By calculating Information Content and Mutual Information over many generations
        of the population, the preferences of Transcription Factor Binding Sites towards
        one manner of storing information can be explored
        
        This method measures MI by iterating through the TFBS list 
        and measuring the frequency of 1-mer and 2-mers for each pair of columns. 
        MI is then computed as the difference between the apriori entropy and the aposteriori entropy:
        
        Let A and B be two random variables or two columns of the motif
        MI(A;B) = H(A) + H(B) - H(A,B)
        
        If MI_mode = 2 then the following is returned. Let A,B,C,D,E,and F be columns of the motif:
        MI(A;B) + MI(C;D) + MI(E;F)
        
        If MI_mode = 1 then the following is returned:
        MI(A;B) + MI(B;C) + MI(C;D) + MI(D;E) + MI(E;F)        
                
        '''
        #Initializes variables
        dinuc_frequency_dic = {}
        firstnuc_frequency_dic = {}
        secondnuc_frequency_dic = {}
        Hx = 0
        Hy = 0
        Hxy = 0
        
        #Sets MI to zero
        self.MI = 0
        
        #Populates dictionaries based on permutations of nucleotides
        for element_outer in self._alphabet:
            firstnuc_frequency_dic[element_outer] = 0
            secondnuc_frequency_dic[element_outer] = 0
            
            for element_inner in self._alphabet:
                dinuc_frequency_dic[element_inner + element_outer] = 0
                
        #Iterates through adjacent pairs with distance self.MI_mode
        for i in range(0,self.length_bs-1,self.MI_mode):
            
            #Fills frequency arrays with dinucleotide and singlenucleotide counts
            for binding_site in self.TFBS:
                dinuc_frequency_dic[binding_site[i] + binding_site[i+1]] += 1
                firstnuc_frequency_dic[binding_site[i]] += 1
                secondnuc_frequency_dic[binding_site[i+1]] += 1
                
            #Converts single nucleotide counts to frequencies
            for key in self._alphabet:
                
                #Converts counts to frequencies
                firstnuc_frequency_dic[key] /= self.number_bs
                secondnuc_frequency_dic[key] /= self.number_bs
                
                #Calculates Hx and Hy
                if firstnuc_frequency_dic[key] != 0:
                    Hx += -(firstnuc_frequency_dic[key] * np.log2(firstnuc_frequency_dic[key]))
                
                if secondnuc_frequency_dic[key] != 0:
                    Hy += -(secondnuc_frequency_dic[key] * np.log2(secondnuc_frequency_dic[key]))
                
                #Sets dictionary to zero
                firstnuc_frequency_dic[key] = 0
                secondnuc_frequency_dic[key] = 0
                
            #Converts dinucleotide counts to frequencies 
            for key in dinuc_frequency_dic.keys():
                
                bp_prob = dinuc_frequency_dic[key] / self.number_bs
                
                #Disregards if frequency is zero
                if bp_prob != 0:
                    
                    #Calculates joint entropy            
                    Hxy += (-1 * (bp_prob * np.log2(bp_prob)))
            
                #Sets the count to zero for the next loop
                dinuc_frequency_dic[key] = 0
                
            #Calculates MI of a pair of columns and adds it to class variable
            self.MI += Hx + Hy - Hxy
                        
            Hx = 0
            Hy = 0
            Hxy = 0
    
    # def Fast_MI(self):
    #     '''
    #     Version of Compute_MI that uses variables generated from Compute_TI and Compute_IC
    #     to calculate MI. Saves computational time.
        
    #     This method takes advantage of the prior computation performed by 
    #     Compute_TI to calculate joint entropy across all adjacent pairs (no overlap)
    #     and Compute_IC to calculate uncertainty of each column. This method then calculates
    #     MI as the difference between cumulative uncertainty and cumulative joint entropy
        
    #     Only usable if MI_mode = 2.
        
    #     DEPRECATED
    #     '''
    #     #Sets MI to zero
    #     self.MI = 0
        
    #     #Converts IC back to uncertainty
    #     cum_entropy = -(self.IC) + 2 * self.length_bs
        
    #     #Uses TI to calculate joint entropy
    #     cum_joint_entropy = 2*(self.length_bs) - self.TI
        
    #     #Calculates MI
    #     self.MI = cum_entropy - cum_joint_entropy
    def Compute_Information(self):
        """
        Organism method used to compute total information stored in an organism's
        list of TFBS. Total Information sums both Information Content and Mutual Information 
        of two columns. In this algorithm, fitting for Total Information allows us to explore 
        how TFBS prefer one manner of storing information over another.
        
        This method is called to update the motif's total information attribute 
        after initialization or after mutation. This method simply computes
        total information as a sum of Information Content and Mutual Information.
        In order to update Information Content and Mutual Information, 
        other methods are called
        """
        
        self.set_IC()
        self.set_MI()
        self.TI = self.IC + self.MI
        
    def Mutate(self):
        '''
        Organism method used to apply mutations at random positions of the organism's 
        genome (The list of TFBS). In order to investigate information conveyed
        column-wise and pair-wise, the ability of information (both MI and IC)
        to withstand mutation must also be explored. This method provides the means
        for mutation to permit this exploration and introduces genetic diversity within 
        the population so as to permit these populations to scale the fitness landscape.
        
        The method uses the mutation rate and genome length to calculate an expected 
        number of mutations, referencing the poisson distribution to get a number 
        of point mutations occuring at this round of mutation for an organism. 
        The number of mutations is used to initialize a list of length (num_of_mut) containing unique 
        positions within the genome for mutation. Iterating through the list,
        the base pair at each position is mutated using a weight based on GC content.
        
        Note: Currently, mutations that mutate a base pair to the same base pair are enabled (Ex: G --> G)
        '''
        
        #Calculates expected value as mt rate * length of sequence * number of mutations
        expected_mut = self.mut_rate * self.number_bs * self.length_bs
        
        #Generates the # of mutations based on expected value using the poisson distribution
        num_mut = np.random.poisson(expected_mut)
        
        #Indicates cutoff value for each base pair
        g_freq = self.GC/2
        c_freq = self.GC
        a_freq = self.GC + 0.5 * (1-self.GC)
        
        #Prepares to draw num_mut amount of numbers between 1 and the genome length (sites for mutation)
        collection_of_positions = []
        sufficient_mut = 0
        
        #Iterative loop to ensure num_mut of unique mutations 
        while sufficient_mut < num_mut: 
            temp_pos = random.randrange(0,(self.number_bs*self.length_bs) - 1)
            
            '''
            #Note: Look for function that sorts and checks progressive elements
            if not(temp_pos in collection_of_positions):
                collection_of_positions.append(temp_pos)
                sufficient_mut += 1
            '''
            sufficient_mut += 1
            collection_of_positions.append(temp_pos)
        
        #Iterates through list and finds referenced position within TFBS
        for pos in collection_of_positions:
            
            #Locates position within list of TFBS
            row = pos // self.length_bs
            col = pos % self.length_bs
            
            #Generates random float between 0 and 1 used to mutate the positions stored in the above array
            #Mutations weighted by GC content
            rand_float = random.random()
            
            #Applies mutation
            if rand_float < g_freq:
                self.TFBS[row] = self.TFBS[row][:col] + 'G' + self.TFBS[row][(col+1):]
            
            elif rand_float < c_freq:
                self.TFBS[row] = self.TFBS[row][:col] + 'C' + self.TFBS[row][(col+1):]
            
            elif rand_float < a_freq:
                self.TFBS[row] = self.TFBS[row][:col] + 'A' + self.TFBS[row][(col+1):]
            
            else:
                self.TFBS[row] = self.TFBS[row][:col] + 'T' + self.TFBS[row][(col+1):]
        
        # Update information metrics
        self.set_IC()
        self.set_MI()
        self.TI = self.IC + self.MI
    
    def Export_Fasta(self,file_name):
        '''
        After the object's constructor is called, the following method can be used
        to export the array of TFBS into a fasta file
        
        Converts array of TFBS into an array of Seq Objects and writes Seq Objects into fasta file. 
        
        Returns
        -------

        '''
        
        #Initializes Seq Object
        seq_list = []
        
        #Creates list of seq objects with sequences corresponding to the organism's TFBS
        for bs in self.TFBS:
            seq_list.append(SeqRecord(Seq(bs),id = ""))
            
        #Writes fasta file if file_name is valid fasta file name.
        if ((file_name[-3:] == 'fas') or (file_name[-5:] == 'fasta')):
            SeqIO.write(seq_list,file_name,'fasta')
        
        else:
            raise(Exception("Extension is atypical for fasta files. Please enter a valid file name with a .fas or .fasta extension name"))
    
    def Generate_Logo(self,file_name):
        '''
        Method to generate logo of organism's TFBS 
        
        Parameters
        -----------
        file_name - Str - Name of the file to be output to (Must be png)
        
        Returns
        -------
        None
        
        '''
        #Initializes Seq Object
        seq_list = []
        
        #Creates list of seq objects with sequences corresponding to the organism's TFBS
        for bs in self.TFBS:
            seq_list.append(Seq(bs))
            
        #Creates motif and generates logo
        object_motif = motifs.create(seq_list)
        
        object_motif.weblogo(file_name, format = 'PNG')
    
    def generate_diPWM(self):
        return create_dipwm_object(self.TFBS)
    
    def all_kmer_scores(self):
        ''' Returns a dictionary where every kmer is associated to its diPWM score. '''
        diPWM = self.generate_diPWM()
        all_kmers = [''.join(kmer) for kmer in itertools.product('ACGT', repeat=self.length_bs)]
        all_kmer_scores = [diPWM.set_score_from_str(kmer) for kmer in all_kmers]
        return dict(zip(all_kmers, all_kmer_scores))
    
        
        



if __name__ == "__main__":
    motifOrg = MotifModel(64, 4,GC_cont = 0.50, mut_rate= 0.01,ID = 1) 
    print(motifOrg,"Motif organism")
    print("Total Information:", motifOrg.TI)
    print("Information Content:",motifOrg.IC)
    print("Mutual Information:",motifOrg.MI)







