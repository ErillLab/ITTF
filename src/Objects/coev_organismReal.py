# -*- coding: utf-8 -*-
"""
Created on Tue May 28 17:05:53 2024

@author: thele
"""

from organism import Organism
from coev_genome import coevGenome, create_dipwm_object, scoreGenome
import copy
import random
from sklearn import metrics
import numpy as np
from statistics import mean

class coevOrganism():
    """
    Class that defines a "coevOrganism" object. The coevolutionary simulation 
    involves evolution through mutation of both the genome and the motif. The 
    objective of the coevolutionary simulation is to gain insight into the 
    fitness landscape of MI and IC motifs in genomes. We hope to prove that
    motifs that encode nearly all of their information through MI with maximum 
    amount of fitness (1 AUPRC) evolve to IC motif. 
    
    This class defines an organism object, which is a genome-motif pair. 
    The class involves methods to mutate the genome and subsequent motif. It can then
    return the subsequent fitness (AUPRC of true binding sites) and information of the motif.
    This object is initialized for each motif in the coevPopulation object
    
    In the coevolutionary simulation, the motif and genome are considered an organism 

    """
    def __init__(self, genomeSize, motif_numberBS, lengthBS, numBS):
        #Initializes attributes 
        self.genomeObj = coevGenome(genomeSize, motif_numberBS, lengthBS, numBS)
        
        #Initializes motif organism from TFBS of genome
        self.BSobject = Organism(numBS, lengthBS,0.5,0,"hi", MI_mode=1, read_in = self.genomeObj.returnBindingSites())
        

    def induceMutation(self):
        """
        Mutations are induced by genome object. Motif (organism object) must be updated by mutations

        Returns
        -------
        None.

        """
        #Induces mutation in genome
        self.genomeObj.induceMutation()
        
        #Updates binding site object and recomputes information
        self.BSobject.TFBS = self.genomeObj.returnBindingSites()
        self.BSobject.Compute_Information()
    
    def returnAUPRC(self):
        """
        Returns area under the precision recall curve of TF (organism) as classifier of genome BS and nonBS

        Returns
        -------
        None.

        """
        #Prepares kulkovsky model for scoring sites using total information of motif
        diP = create_dipwm_object(self.BSobject.TFBS)
        
        #Using Kulakovsky model, scores circular genome using TF made from organism's motif
        scoresList =  scoreGenome(diP, self.genomeObj.lengthBS, self.genomeObj.getCircularGenome())
        
        #Creates a list the length of the genome with 1 for positive (BS) and 0 for negatives (everything else)
        PNlist = []
        positiveIndices = self.genomeObj.indicesBS
        for i in range(self.genomeObj.genomeSize):
            
            #Adds 1 if index is a positive (binding site). Otherwise adds 0
            if (i in positiveIndices):
                PNlist.append(1)
            else:
                PNlist.append(0)
        
        #Computes AUC
        precisionArray, recallArray, thresholds = metrics.precision_recall_curve(PNlist, scoresList)
        organismAUC = metrics.auc(recallArray, precisionArray)
        
        return organismAUC
    
    def organismInformation(self):
        """
        Returns Information (first element, second element, third element) = 
        (IC, MI, TI)

        Returns
        -------
        None.

        """
        return (self.BSobject.IC, self.BSobject.MI, self.BSobject.TI)
        
        
        
class coevPopulation():
    """
    Population composed of coevOrganisms. Used to run the coevolutionary simulation
    to coevolve binding sites with motifs. This population is both 
    used to generate MI/IC coevOrganisms and run the experiments to determine
    the mutational stability of MI coevOrganisms
    
    The fitness function for evolving the population is a linear combination 
    of the AUPRC (TF's classification ability) and %MI deviation 
                  
    Each generation, the following cycle occurs:
        Replication with # of offspring based on fitness -> Mutation of all organisms -> Random sampling of popSize # of organisms 
    popSize (int) - Size of population (# of coevOrganisms)
    targetMIproportion (float) - Goal proportion MI for population. Involved in the fitness function for evolving the population
    Note: Should be bounded from 0 to 1
    
    alpha (int) - Tunable parameter weighting for target_MIpercent (alpha = 1 means target_MIpercent is weighted same as AUPRC)
    beta (int) - Weight for the number of offspring for each organism. 
    idList (list) - List of strings used to keep track of which offspring come from which ancestors in the original population.
    
    """
    def __init__(self, popSize, targetMIproportion, alpha, beta, genomeSize, motif_numberBS, lengthBS, numBS, idList = None):
        #Initializes attributes
        self.popSize = popSize
        self.targetMIproportion = targetMIproportion
        self.alpha = alpha
        self.beta = beta
        
        #Throws error if targetMIproportion is not bounded from 0 to 1
        if (self.targetMIproportion >= 1) | (self.targetMIproportion <= 0):
            raise ValueError("targetMIproportion is a proportion not a percent! Bound from 0 to 1")
        
        #Defines pertient lists
        self.populationList = [] #Keeps track of all organisms
        self.fitnessList = [] #Keeps track of the fitness of all organisms
        self.replicatedList = [] #Keeps track of all organisms after replication 
        
        #Initializes population
        for i in range(self.popSize):
            self.populationList.append(coevOrganism(genomeSize, motif_numberBS, lengthBS, numBS))
        
        #Ensures idList is either none or the size of the populationList
        if ((idList != None) & (len(idList) != self.popSize)):
            raise ValueError("IDlist is not same size as population List in coevolutionary algorithm")
        else:
            #Stores ID in ancestor organism
            for i in range(self.popSize):
                self.populationList[i].orgID = idList[i]
    
    def determineFitness(self):
        self.fitnessList = []
        
        #Initializes list for keeping track of % information
        IClist = []
        MIlist = []
        TIlist = []
        AUPRClist = []
        
        self.maxOrganismFitness = 0
        
        #Iterates through list and determines fitness of each organism
        for org in self.populationList:      
            #Grabs AUPRC
            orgAUPRC = org.returnAUPRC()
            
            #Grabs %MI
            informationTuple = org.organismInformation()
            orgProportionMI = informationTuple[1]/informationTuple[2] #MI/TI
            
            #Defines fitness function
            MIdeviation = abs(self.targetMIproportion - orgProportionMI)
            
            #Computes org fitness and adds to fitness list. Bounded from 0 to 1
            orgFitness = ((self.alpha * (1 - MIdeviation)) + orgAUPRC)/(1+self.alpha)
            self.fitnessList.append(orgFitness)
            
            #Keeps track of maximum Fitness in population
            if (orgFitness > self.maxOrganismFitness):
                self.maxOrganismFitness = orgFitness
                 
            #Keeps track of information
            IClist.append(informationTuple[0])
            MIlist.append(informationTuple[1])
            TIlist.append(informationTuple[2])
            AUPRClist.append(orgAUPRC)
            
        #Variable for keeping track of mean fitness
        self.meanFitness = mean(self.fitnessList)
        
        #Variable for keeping track of mean Information
        self.meanInformation = (mean(IClist), mean(MIlist), mean(TIlist))
        
        self.meanAUPRC = mean(AUPRClist)
    
    def Replicate(self):
        """
        Ran after determineFitness. Computes # of organisms based on fitness and replicates

        Returns
        -------
        None.

        """
        self.replicatedList = []
        
        #Defines mean for poisson distribution. Computed as exponentiated difference from mean fitness
        r = lambda s: self.beta **(s-self.meanFitness)  
        
        #Iterates through population, determining # of offspring
        for i in range(len(self.populationList)):
            
            #Determines mean offspring as a function of fitness for current offspring
            meanOffspring = r(self.fitnessList[i])
            numOffspring = np.random.poisson(meanOffspring)
            
            #Adds current parent organism to replicatedList
            self.replicatedList.append(self.populationList[i])
            
            #Adds offspring to replicatedList
            for j in range(numOffspring):
                self.replicatedList.append(copy.deepcopy(self.populationList[i]))
        
            
    def Mutate(self):
        """
        Mutates entire population after reproduction. Note: We must call this method after Replicate

        Returns
        -------
        None.

        """
        for i in (range(len(self.replicatedList))):
            self.replicatedList[i].induceMutation()
    
    def Select(self):
        """
        Selection of popSize # of organisms from replicatedList. Must be called after Mutate method
        """
        #Reinitializes populatonList and fitnessList due to selection
        self.populationList = []
        self.fitnessList = []
        
        #First raises error if replicatedList is ever smaller than popSize
        if (len(self.replicatedList) < self.popSize):
            raise ValueError("Replicated List is smaller than popSize")
            
        #Shuffles replicatedList and puts first 100 randomly selected organisms in populationList
        random.shuffle(self.replicatedList)
        self.populationList = self.replicatedList[:self.popSize]
        
        #Reinitializes replicatedList
        self.replicatedList = []
    
    def runOneGeneration(self):
        """
        Helper function for coevolutionary simulation that runs 
        one generation of the simulation:
        Determines Fitness, Reproduces, Mutates, and Selects
        
        Note: runOneGeneration is meant to be used in an iterative loop

        Returns
        -------
        None.

        """
        #Begins with a population of size: Pop size
        self.determineFitness()
        self.Replicate()
        self.Mutate()
        self.Select()
        #Ends with a population of size: Pop size. 

                    
            
            
        
        
        
            
        
            
        
            
        
        
        
        
    
        