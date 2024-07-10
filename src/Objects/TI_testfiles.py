# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:58:58 2023

@author: thele
"""
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random

def write_fasta (type_seq,num_seq,len_seq):
    seq_list = []
    alphabet = ['A','C','G','T']
    
    if type_seq == "perfect conservation":
        for i in range(num_seq):
            seq_list.append(SeqRecord(Seq(''.join(['A' for i in range(len_seq)]))))
        
        SeqIO.write(seq_list,'../TI_test/perfect_conservation.fas','fasta')
    
    elif type_seq == "completely random":
        for i in range(num_seq):
            seq_list.append(SeqRecord(Seq(''.join(random.choices(alphabet, k = num_seq)))))
        
        SeqIO.write(seq_list,'../TI_test/completely_random.fas','fasta')
            
    elif type_seq == "max_MI":
        
        sequenceA = "ACGTCGTGAGCG"
        sequenceB = "TGCAGCACTCGC"
        sequenceC = "GATGATGAGATT"
        sequenceD = "CTACTACTCTAA"
        
        for i in range(250):
            seq_list.append(SeqRecord(Seq(sequenceA)))
            seq_list.append(SeqRecord(Seq(sequenceB)))
            seq_list.append(SeqRecord(Seq(sequenceC)))
            seq_list.append(SeqRecord(Seq(sequenceD)))



            
            
        SeqIO.write(seq_list,'../TI_test/max_MI.fas','fasta')
                
write_fasta('max_MI',1000,20)