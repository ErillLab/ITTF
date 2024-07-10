from Bio import SeqIO
import os
import pickle
from math import log2
from statistics import mean


def computeMixtureEntropy(realBS):
    """
    Takes in list of real Binding Sites and computes the KLD of binding site frequency distribution to uniform (Entropy)

    Parameters
    ----------
    realInstances : List of binding sites

    Returns
    -------
    Mixture Entropy

    """
    alphabet = ["A","C","G","T"]
    uniformDistribution = {}
    realDistribution = {}
    
    for alphabet1 in alphabet:
        for alphabet2 in alphabet:
            realDistribution[alphabet1 + alphabet2] = 0
            uniformDistribution[alphabet1 + alphabet2] = 1/16
    
    for bs in realBS:
        realDistribution[bs] += 1
    sumCounts = sum(realDistribution.values())
    for key in realDistribution.keys():
        realDistribution[key] = realDistribution[key] / sumCounts
        
    #Begins computing Mixture Entropy
    MixtureEntropy = 0
    for key in realDistribution.keys():
        MixedProportion = (1/2) * (realDistribution[key] + uniformDistribution[key])
        MixtureEntropy += (MixedProportion * log2(MixedProportion))
    
    return (-1 * MixtureEntropy)
        
randomMotifs = [instances.seq for instances in list(SeqIO.parse("random_noresidual.fas", format = "fasta"))]
populationFile = "./Output/Populations/Divergence_Experiment1/population_100_128_2_0IC.pickle"

with open(populationFile,'rb') as populationFileVar:
    population = pickle.load(populationFileVar)

MixtureEntropyList = []
for organism in population.population_list:
    MixtureEntropyList.append(computeMixtureEntropy(organism.TFBS))
    
print("0% IC Mixture Distribution Entropy:",mean(MixtureEntropyList))

