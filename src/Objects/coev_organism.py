# -*- coding: utf-8 -*-

from src.Objects.organism import Organism
from src.helperFunctions import create_dipwm_object
from src.Objects.coev_genome import coevGenome, scoreGenome
import copy
import random
from sklearn import metrics
import numpy as np
from statistics import mean, stdev

class coevOrganism():
    """
    Class that defines a "coevOrganism" object. The coevolutionary simulation 
    involves evolution through mutation of both the genome and the motif (recognizre).  
    The objective of the coevolutionary simulation is to gain insight of the
    fitness landscape of MI and IC motifs in genomes. 
    
    This class defines an organism object, which is a genome-motif pair. 
    The class involves methods to mutate the genome and subsequent motif. It can then
    return the subsequent fitness (AUPRC of true binding sites) and information of the motif.
    This object is used to create a population of coevOrganisms in coevPopulation
        
    In the coevolutionary simulation, the motif and genome are considered an organism 
    """
    def __init__(self, genomeSize, motif_numberBS, lengthBS, numBS):
        """
        Constructor to initialize the genome, motif object, and the binding 
        site object. The genome is initialized as a random collection of bases,
        the motif object is initialized from the gene at the beginning of the genome,
        and the binding site object is initialized from the binding sites 
        throughout the genome

        Parameters
        ----------
        genomeSize : int
            Length of genome, including motif and binding sites. Must be long 
            enough to include motif length and 2 * length of all binding sites
        motif_numberBS : int
            Number of binding sites in gene. Higher values lead to more
            fine-grained information.
        lengthBS : int
            Length of each binding site.
        numBS : int
            Number of binding sites in genome
        """
        #Initializes attributes 
        self.genomeObj = coevGenome(genomeSize, motif_numberBS, lengthBS, numBS)
        
        #Initializes motif object from gene in genomeObj
        self.geneObj = self.genomeObj.geneObj
        
        #Initializes binding site organism from TFBS of genome
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
        
        #Updates motif gene
        self.geneObj = self.genomeObj.geneObj
        
        #Updates binding site object and recomputes information
        self.BSobject.TFBS = self.genomeObj.returnBindingSites()
        self.BSobject.Compute_Information()
        
    def induceMutationFrozenRecognizer(self):
        """
        Mutations are induced by genome object. Recognizer remains frozen (non mutated)

        Returns
        -------
        None.

        """
        #Induces mutation in genome
        self.genomeObj.induceMutation_freezeRecognizer()
        
        #Updates motif gene
        self.geneObj = self.genomeObj.geneObj
        
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
        diP = create_dipwm_object(self.geneObj.transcribeGene())
        
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
    
    def returnMotifInformation(self):
        """
        Returns Information of recognizer (first element, second element, third element) = 
        (IC, MI, TI)

        Returns
        -------
        Tuple of information: (IC, MI, TI)
        
        """
        return self.geneObj.returnMotifInformation()
    
    def returnBSInformation(self):
        """
        Returns Information of binding sites (first element, second element, third element) = 
        (IC, MI, TI)

        Returns
        -------
        Tuple of information: (IC, MI, TI)

        """
        return (self.BSobject.IC, self.BSobject.MI, self.BSobject.TI)
    
    def returnJSD(self):
        """
        Returns Jennson Shannon divergence between the motif and the binding sites
        
        Jennson Shannon Divergence function computes the divergence between 
        frequency distributions with two base pairs. As such, JSD between the 
        motif and binding site frequency distributions are computed for a pair 
        of columns. A sliding window is used to compute the total JSD by 
        summing JSD of successive columns
        
        ???WIP
        
        Returns
        -------
        totalJSD: float
            sum of the JSD computed between successive columns of the motif 
            and the binding sites
        """
        
        motifSites = self.geneObj.transcribeGene()
        bindingSites = self.BSobject.TFBS
        
        #lengthBS
        

        
        
        
        
        
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
        if (idList != None):
            if ((len(idList) != self.popSize)):
                raise ValueError("IDlist is not same size as population List in coevolutionary algorithm")
            else:
                #Stores ID in ancestor organism
                for i in range(self.popSize):
                    self.populationList[i].orgID = idList[i]
    
    def determineFitness(self):
        self.fitnessList = []
        
        #Iterates through list and determines fitness of each organism
        for org in self.populationList:      
            #Grabs AUPRC
            orgAUPRC = org.returnAUPRC()
            
            #Grabs %MI
            informationTuple = org.returnMotifInformation()
            orgProportionMI = informationTuple[1]/informationTuple[2] #MI/TI
            
            #Defines fitness function
            MIdeviation = abs(self.targetMIproportion - orgProportionMI)
            
            #Computes org fitness and adds to fitness list. Bounded from 0 to 1
            orgFitness = ((self.alpha * (1 - MIdeviation)) + orgAUPRC)/(1+self.alpha)
            self.fitnessList.append(orgFitness)
        
                             
    
    def Replicate(self):
        """
        Ran after determineFitness. Computes # of organisms based on fitness and replicates

        Returns
        -------
        None.

        """
        self.replicatedList = []
        self.meanFitness = mean(self.fitnessList)
        
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
    
    def generationFitness(self):
        """
        Called in runOneGeneration to return information of each organism 
        in the population after a generation of selection

        Returns
        -------
        popDict - Dict
            Dictionary with information about the population. Contains 
            elements describing the mean, standard deviation, and max
            for the following: motif information, AUPRC, and fitness 

        """       
        #Summary variable of population
        popDict = dict()
        
        #Initializes list for keeping track of % information
        IClist = []
        MIlist = []
        TIlist = []
        AUPRClist = []
        
        #Variable for maximumOrganismFitness
        maxOrganismFitness = 0
        fitnessList = []
        
        #Iterates through list and determines fitness of each organism
        for org in self.populationList:      
            #Grabs AUPRC
            orgAUPRC = org.returnAUPRC()
            
            #Grabs %MI
            informationTuple = org.returnMotifInformation()
            orgProportionMI = informationTuple[1]/informationTuple[2] #MI/TI
            
            #Defines fitness function
            MIdeviation = abs(self.targetMIproportion - orgProportionMI)
            
            #Computes org fitness and adds to fitness list. Bounded from 0 to 1
            orgFitness = ((self.alpha * (1 - MIdeviation)) + orgAUPRC)/(1+self.alpha)
            fitnessList.append(orgFitness)
            
            #Keeps track of maximum Fitness in population
            if (orgFitness > maxOrganismFitness):
                maxOrganismFitness = orgFitness
                 
            #Keeps track of information
            IClist.append(informationTuple[0])
            MIlist.append(informationTuple[1])
            TIlist.append(informationTuple[2])
            AUPRClist.append(orgAUPRC)
            
        #Variable for keeping track of mean fitness
        popDict["Mean Fitness"] = mean(fitnessList)
        if self.popSize > 1:
            popDict["Stdev Fitness"] = stdev(fitnessList)
        popDict["Top Fitness"] = max(fitnessList)
        
        #Variable for keeping track of mean Information
        popDict["Mean Information"] = (mean(IClist), mean(MIlist), mean(TIlist))
        if self.popSize > 1:
            popDict["Stdev Information"] = (stdev(IClist), stdev(MIlist), stdev(TIlist))
        
        popDict["Mean AUPRC"] = mean(AUPRClist)
        if self.popSize > 1:
            popDict["Stdev AUPRC"] = stdev(AUPRClist)
        popDict["Top AUPRC"] = max(AUPRClist)
        
        return popDict
        
        
    def runOneGeneration(self):
        """
        Helper function for coevolutionary simulation that runs 
        one generation of the simulation:
        Determines Fitness, Reproduces, Mutates, and Selects
        
        Note: runOneGeneration is meant to be used in an iterative loop

        Returns
        -------
        self.generationFitness() - Returns a dictionary with information 
        about the population

        """
        #Begins with a population of size: Pop size
        self.determineFitness()
        self.Replicate()
        self.Mutate()
        self.Select()
        
        #Returns dictionary with information about the population
        return self.generationFitness()
    
                
    def Mutate_FrozenRecognizer(self):
        """
        Mutates entire population after reproduction. Note: We must call this method after Replicate

        Returns
        -------
        None.

        """
        for i in (range(len(self.replicatedList))):
            self.replicatedList[i].induceMutationFrozenRecognizer()
            
    def runOneGeneration_FrozenRecognizer(self):
        """
        Helper function for experiments to simulate the longer 
        evolutionary timeline of the motif. Carries out one generation of 
        the evolutionary simulation without mutating the motif gene (frozen recognizer)
    

        Returns
        -------
        None.

        """
        #Begins with a population of size: Pop size
        self.determineFitness()
        self.Replicate()
        self.Mutate_FrozenRecognizer()
        self.Select()
        
        #Returns dictionary with information about the population
        return self.generationFitness()
    

                    
            
            
if __name__ == "__main__":
    """
    Generates a sample coevolutionary population
    """
    exampleCoevPop = coevPopulation(1,0.50,1,10,2046,32,8,32)
    print(exampleCoevPop.runOneGeneration())
    
        
        
            
        
            
        
            
        
        
        
        
    
        