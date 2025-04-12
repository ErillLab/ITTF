# -*- coding: utf-8 -*-

from src.Objects.motiforganism import MotifModel
from src.helperFunctions import create_dipwm_object
from src.Objects.genome import Genome, scoreGenome
import copy
import random
from sklearn import metrics
import numpy as np
from statistics import mean, stdev
import json
import sys
import heapq
import operator

class coevOrganism():
    """
    Class that defines a "coevOrganism" object. The coevolutionary simulation 
    involves evolution through mutation of both the genome and the motif (recognizre).  
    The objective of the coevolutionary simulation is to gain insight of the
    fitness landscape of MI and IC motifs in genomes. 
    
    This class defines an organism object, which is a genome-recognizer pair. 
    The class involves methods to mutate the genome and subsequent recognizer. It can then
    return the subsequent fitness (AUPRC of true binding sites) and information of the recognizer.
    This object is used to create a population of coevOrganisms in coevPopulation
        
    In the coevolutionary simulation, the motif and genome are considered an organism 
    """
    def __init__(self, genomeSize = None, motif_numberBS = None, lengthBS = None, numBS = None, idTag = None,importJSON = None):
        """
        Constructor to initialize the genome, recognizer object, and the binding 
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
        idTag: str
            Optional tag used to track ancestors of current organism in coevolutionary simulation
        importJSON: str/Dict/None
            Used to initialize object from json file path, python dict object.
            Set to None to initialize normally. Other variables become unnecessary to specify
        """
        #Initialze Normally
        if importJSON is None:
            #Initializes attributes 
            self.genomeObj = Genome(genomeSize, motif_numberBS, lengthBS, numBS)
            
            #Initializes recognizer object from gene in genomeObj
            self.recognizerObj = MotifModel(motif_numberBS, lengthBS,0.5,0,"Deprecated", MI_mode=1, read_in = self.genomeObj.transcribeGene())
            
            #Initializes binding site organism from TFBS of genome
            self.BSobject = MotifModel(numBS, lengthBS, 0.5,0,"Deprecated", MI_mode=1, read_in = self.genomeObj.returnBindingSites())
            
            self.orgID = idTag
            
        #Initialze from file or object
        else:
            #From file
            if type(importJSON) is str:
                self._import(filepath = importJSON)
            
            #From object
            else:
                self._import(importObject = importJSON)
        

    def induceMutation(self,mutateGenome = True,frozenRecognizer = False):
        """
        Mutations are induced by genome object. 
        
        If frozenRecognizer is set to false, whole genome is mutated. Otherwise,
        recognizer is "frozen" and mutation hits the rest of the genome

        Returns
        -------
        None.

        """
        #Induces mutation in genome but freezes recognizer
        self.genomeObj.induceMutation(mutateGenome, frozenRecognizer)
        
        #Updates recognizer (only necessary if frozenRecognizer = False)
        self.recognizerObj.TFBS = self.genomeObj.transcribeGene()
        self.recognizerObj.Compute_Information()
        
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
        diP = create_dipwm_object(self.recognizerObj.TFBS)
        
        #Using Kulakovsky model, scores circular genome using TF made from organism's motif
        scoresList =  scoreGenome(diP, self.genomeObj.lengthBS, self.genomeObj.getCircularGenome())
        
        positiveIndices = self.genomeObj.indicesBS

        #Creates a list the length of the genome with 1 for positive (BS) and 0 for negatives (everything else)
        PNlist = np.zeros(self.genomeObj.genomeSize).astype(int)
        PNlist[positiveIndices] = 1
        
        #Computes AUC
        precisionArray, recallArray, thresholds = metrics.precision_recall_curve(PNlist, scoresList)       
        return metrics.auc(recallArray, precisionArray)
    
    def returnRecognizerInformation(self):
        """
        Returns Information of recognizer (first element, second element, third element) = 
        (IC, MI, TI)

        Returns
        -------
        Tuple of information: (IC, MI, TI)
        
        """
        return (self.recognizerObj.IC, self.recognizerObj.MI, self.recognizerObj.TI)
    
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
        
    
    def _import(self, filepath = None, importObject = None):
        """
        Imports object from json file at filepath or from python object

        Parameters
        ----------
        filepath : str or None
            Path to json file for importing. If none, must specify python object in importObject
        importObject: Dict or None
            By default, none if filepath is specified. If filepath is not specified, should be set to python Dict
        Returns
        -------
        None.

        """
        #Initializes genome obj from file 
        if filepath:
            self.genomeObj = Genome(importJSON = filepath)
            
            #Opens file to grab ID (coevOrganism unique attribute)
            with open(filepath, 'r') as f:
                temp_dic = json.load(f)
                self.orgID = temp_dic["orgID"]
                f.close()
            
        #Otherwise initializes genome obj from object
        else:
            self.genomeObj = Genome(importJSON = importObject)
            
            #Grabs orgID (coevOrg unique attribute)
            self.orgID = importObject['orgID']
            
        
        #Grabs attributes of genomeObj
        import_dic = self.genomeObj._export()
        
        #Initializes attributes 
        #??? Rewrite when import/export functions are made for Motif Model/Population classes
        self.recognizerObj = MotifModel(N = import_dic['motif_numberBS'], 
                                        L = import_dic['lengthBS'], 
                                        GC_cont = 0.5, 
                                        mut_rate = 0, 
                                        ID = "Deprecated",
                                        MI_mode = 1,
                                        read_in = self.genomeObj.transcribeGene())
        
        self.BSobject = MotifModel(N = import_dic['numberBS'], 
                                   L = import_dic['lengthBS'], 
                                   GC_cont = 0.5, 
                                   mut_rate = 0, 
                                   ID = "Deprecated",
                                   MI_mode = 1, 
                                   read_in = self.genomeObj.returnBindingSites())
        
        
    def _export(self, filepath = None):
        """
        Exports object to json file at filepath or to python Object

        Parameters
        ----------
        filepath : str or None
            Path to json file for exporting. If none, returns python object
            
        Returns
        -------
        None or exported object
        """
        
        #Calls export function and adds additional attribute
        out_dict = self.genomeObj._export()
        out_dict['orgID'] = self.orgID
        
        #Prints to file or returns
        if filepath is not None:
            with open(filepath, 'w') as f:
                json.dump(out_dict, f)
                
        #Otherwise returns object
        else:
            return out_dict
        
            
        

        
        
        
        
        
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
    targetMIproportion (float) - Goal proportion MI for population. Involved in the fitness function for evolving the population. 
    Note: Should be bounded from 0 to 1
    targetAUPRC (float) - Goal AUPRC for population. Involved in the fitness function for evolving the population. 
    Note: Should be bounded from 0 to 1
    alpha (int) - Tunable parameter weighting for target_MIpercent (alpha = 1 means target_MIpercent is weighted same as AUPRC)
    beta (int) - Weight for the number of offspring for each organism. 
    idList (list) - List of strings used to keep track of which offspring come from which ancestors in the original population.
    importList (list or None) - List of coevOrganisms used to initialize populationList. Organisms are duplicated to fill populationList
    Note: None will randomly initialize population list
    importFile (str or None) - JSON file path to import from. If specified, other inputs are unnecessary  
    """
    def __init__(self, popSize = None, targetMIproportion = None, targetAUPRC = None, alpha = None, beta = None, genomeSize = None, motif_numberBS = None, lengthBS = None, numBS = None, idList = None, importList = None, importFile = None):
        
        #Initialize normally
        if importFile is None:
            #Initializes attributes
            self.popSize = popSize
            self.targetMIproportion = targetMIproportion
            self.targetAUPRC = targetAUPRC
            self.alpha = alpha
            self.beta = beta
            
            #Throws error if targetMIproportion is not bounded from 0 to 1
            if (self.targetMIproportion > 1) | (self.targetMIproportion < 0):
                raise ValueError("targetMIproportion is a proportion not a percent! Bound from 0 to 1")
            
            #Defines pertient lists
            self.populationList = [] #Keeps track of all organisms
            self.fitnessList = [] #Keeps track of the fitness of all organisms
            self.replicatedList = [] #Keeps track of all organisms after replication 
            

            #Initializes population 
            for i in range(self.popSize):
                
                #Initializes population with passed organisms
                if importList:
                    self.populationList.append(copy.deepcopy(importList[i % len(importList)]))
                    
                else:
                    #Otherwise randomly initializes organisms
                    self.populationList.append(coevOrganism(genomeSize, motif_numberBS, lengthBS, numBS))
            
            #Ensures idList is either none or the size of the populationList
            if (idList != None):
                if ((len(idList) != self.popSize)):
                    raise ValueError("IDlist is not same size as population List in coevolutionary algorithm")
                else:
                    #Stores ID in ancestor organism
                    for i in range(self.popSize):
                        self.populationList[i].orgID = idList[i]
        
        #Initialize from file
        else:
            self._import(filepath = importFile)
            
    def returnError(self,orgAUPRC,orgMI):
        """
        Helper function for determineFitness
        
        Computes the error of a passed organism's AUPRC/MI based on target
        as a weighted L2 norm
        
        Returns
        -------
        Weighted eucledian distance (tuple)

        """
        return ((orgAUPRC-self.targetAUPRC)**2 + (self.alpha*(orgMI-self.targetMIproportion))**2)**(1/2)
    
    def determineFitness(self):
        """
        Computes fitness of each organism and saves to self.fitnessList
        
        
        Fitness Function: Fitness = alpha * (Deviation from desired MI) + AUPRC

        Returns summary statistics of population organisms

        Returns
        -------
        popDict - Dict
            Dictionary with information about the population. Contains 
            elements describing the mean, standard deviation, and max
            for the following: motif information, AUPRC, and fitness 
        """
        self.fitnessList = []
        
        #Summary variable of population
        popDict = dict()
        
        #Initializes list for keeping track of % information
        IClist_recognizer = []
        MIlist_recognizer = []
        TIlist_recognizer = []
        percentMI_recognizer = []
        
        IClist_BS = []
        MIlist_BS = []
        TIlist_BS = []
        percentMI_BS = []
        
        AUPRClist = []
        
        #Adjusts ids
        for i in range(self.popSize):
            self.populationList[i].orgID = "org " + str(i)
            
        #Iterates through list and determines fitness of each organism
        for org in self.populationList:      
            #Grabs AUPRC
            orgAUPRC = org.returnAUPRC()
            
            #Grabs %MI
            informationTuple_bs = org.returnBSInformation()
            orgProportionMI_bs = informationTuple_bs[1]/informationTuple_bs[2] #MI/TI
            
            #Keeps track of binding site information
            informationTuple_recognizer = org.returnRecognizerInformation()
            
            
            #Computes org fitness as negative error and adds to fitness list.
            orgFitness = -self.returnError(orgAUPRC = orgAUPRC, orgMI = orgProportionMI_bs)
            self.fitnessList.append(orgFitness)

            #Keeps track of information
            IClist_recognizer.append(informationTuple_recognizer[0])
            MIlist_recognizer.append(informationTuple_recognizer[1])
            TIlist_recognizer.append(informationTuple_recognizer[2])
            percentMI_recognizer.append(informationTuple_recognizer[1]/informationTuple_recognizer[2])
            
            IClist_BS.append(informationTuple_bs[0])
            MIlist_BS.append(informationTuple_bs[1])
            TIlist_BS.append(informationTuple_bs[2])
            percentMI_BS.append(informationTuple_bs[1]/informationTuple_bs[2])
            
            AUPRClist.append(orgAUPRC)
        
        ####Saves descriptive lists to dict
        popDict["Fitness List"] = self.fitnessList
        popDict["AUPRC List"] = AUPRClist
        popDict["Recognizer IC"] = IClist_recognizer
        popDict["Recognizer MI"] = MIlist_recognizer
        popDict["Recognizer TI"] = TIlist_recognizer
        popDict["Recognizer %MI"] = percentMI_recognizer
        popDict["BS IC"] = IClist_BS
        popDict["BS MI"] = MIlist_BS
        popDict["BS TI"] = TIlist_BS
        popDict["BS %MI"] = percentMI_BS
        
        return popDict
        
        

        """
        popDict["Mean Fitness"] = mean(self.fitnessList)
        if self.popSize > 1: 
            popDict["Stdev Fitness"] = stdev(self.fitnessList)
        popDict["Top Fitness"] = max(self.fitnessList)
        
        #Variable for keeping track of mean Recognizer Information
        popDict["Mean Recognizer IC"] = mean(IClist_recognizer)
        popDict["Mean Recognizer MI"] = mean(MIlist_recognizer)
        popDict["Mean Recognizer TI"] = mean(TIlist_recognizer)
        popDict["Mean Recognizer %MI"] = mean(percentMI_recognizer) 

        if self.popSize > 1:
            popDict["Stdev Recognizer IC"] = stdev(IClist_recognizer)
            popDict["Stdev Recognizer MI"] = stdev(MIlist_recognizer)
            popDict["Stdev Recognizer TI"] = stdev(TIlist_recognizer)
            popDict["Stdev Recognizer %MI"] = stdev(percentMI_recognizer)
        
        
        #Computes % Recognizer 
        #Variable for keeping track of mean binding site Information
        popDict["Mean BS IC"] = mean(IClist_BS)
        popDict["Mean BS MI"] = mean(MIlist_BS)
        popDict["Mean BS TI"] = mean(TIlist_BS)
        popDict["Mean BS %MI"] = mean(percentMI_BS)
        if self.popSize > 1:
            popDict["Stdev BS IC"] = stdev(IClist_BS)
            popDict["Stdev BS MI"] = stdev(MIlist_BS)
            popDict["Stdev BS TI"] = stdev(TIlist_BS)
            popDict["Stdev BS %MI"] = stdev(percentMI_BS)
        
            
        popDict["Mean AUPRC"] = mean(AUPRClist)
        if self.popSize > 1:
            popDict["Stdev AUPRC"] = stdev(AUPRClist)
        popDict["Top AUPRC"] = max(AUPRClist)
        
        #Adjusts
        popDict["Parent organisms"] = [org.orgID for org in self.populationList]
        
        return popDict
    """
                             
    
    def Replicate(self):
        """
        Ran after determineFitness. Computes # of organisms based on fitness and replicates
        
        numOffspring = beta^(Fitness of organism - Mean Fitness of Population)
        Returns
        -------
        None.

        """
        sys.stdout.flush()
        self.replicatedList = []
        self.meanFitness = mean(self.fitnessList)
        
        #Defines mean for poisson distribution. Computed as exponentiated difference from mean fitness
        r = lambda s: self.beta **(s-self.meanFitness)  
        
        #Computes ranks of organisms based on fitness
        #First index of tuple is Rank and second is index
        #rankIndex = list(enumerate(list(zip(*sorted(enumerate(self.fitnessList), 
        #                                            key=operator.itemgetter(1), 
        #                                            reverse=True)))[0],start=1))
        
        #Generates list of ranks with same order as self.populationList
        #rankList = [sorted(rankIndex, key = operator.itemgetter(1))[i][0] for i in range(len(rankIndex))]

        #Iterates through population, determining # of offspring
        for i in range(len(self.populationList)):
            
            #Determines mean offspring as a function of rank
            meanOffspring = r(self.fitnessList[i])
            #meanOffspring = self.beta * (1/rankList[i])
            numOffspring = np.random.poisson(meanOffspring)
            
            #Adds current parent organism to replicatedList
            self.replicatedList.append(self.populationList[i])
            
            #Adds offspring to replicatedList
            for j in range(numOffspring):
                self.replicatedList.append(copy.deepcopy(self.populationList[i]))
        
        sys.stdout.flush()
            
    def Mutate(self, mutateGenome = True,frozenRecognizer = False):
        """
        Mutates entire population after reproduction. Note: We must call this method after Replicate
        
        Parameters:
            frozenRecognizer: boolean
                True to mutate the entire genome or False to mutate everything other than the recognizer
                
        Returns
        -------
        None.

        """
        sys.stdout.flush()
        #Mutates every replicated organism
        for i in (range(len(self.replicatedList))):
            self.replicatedList[i].induceMutation(mutateGenome,frozenRecognizer)
        
        sys.stdout.flush()
    
    def Select(self):
        """
        Selection of popSize # of organisms from replicatedList. Must be called after Mutate method
        """
        sys.stdout.flush()
        #Reinitializes populatonList due to selection
        self.populationList = []
        
        
        #First raises error if replicatedList is ever smaller than popSize
        if (len(self.replicatedList) < self.popSize):
            raise ValueError("Replicated List is smaller than popSize")
            
        #Selects most fit organisms for next population  
        random.shuffle(self.replicatedList)
        self.populationList = self.replicatedList[:self.popSize]
        #mostFitOrganismsIndx = list(zip(*heapq.nlargest(self.popSize, 
        #                    enumerate(self.fitnessList), key=operator.itemgetter(1))))[0]
        #self.populationList = [self.replicatedList[i] for i in mostFitOrganismsIndx]
            
        #Reinitializes replicatedList and fitnessList
        self.replicatedList = []
        self.fitnessList = []
        
        sys.stdout.flush()
   
        
    def runOneGeneration(self, mutateGenome = False,frozenRecognizer = False):
        """
        Helper function for coevolutionary simulation that runs 
        one generation of the simulation:
        Determines Fitness, Reproduces, Mutates, and Selects
        
        Used to fit populations to generate 1AUPRC organism and to 
        run the coevolutionary competition experiment
        
        Note: runOneGeneration is meant to be used in an iterative loop
        
        Parameters:
            frozenRecognizer: boolean
                True to mutate the entire genome or False to mutate everything other than the recognizer
                
        Returns
        -------
        self.generationFitness() - Returns a dictionary with information 
        about the population

        """
        #Begins with a population of size: Pop size
        popDict = self.determineFitness()
        self.Replicate()
        
        #Mutates population
        self.Mutate(mutateGenome=mutateGenome,frozenRecognizer=frozenRecognizer)
            
        self.Select()
        
        #Returns dictionary with information about the population
        return popDict
    
    def _import(self,filepath):
        """
        Imports object from json file

        Parameters
        ----------
        filepath : str or None
            Path to json file for importing. 
            
        Returns
        -------
        None.

        """
        with open(filepath, 'r') as f:
            import_dic = json.load(f)
            
        #Initiailizes attributes
        self.popSize = import_dic['populationSize']
        self.targetMIproportion = import_dic['targetMIproportion']
        self.targetAUPRC = import_dic['targetAUPRC']
        self.alpha = import_dic['alpha']
        self.beta = import_dic['beta']
        
        #Initializes populationlist
        self.populationList = [coevOrganism(importJSON = import_dic['populationList'][i])
                                            for i in range(self.popSize)]
        
        
    def _export(self,filepath):
        """
        Exports object to json file at filepath

        Parameters
        ----------
        filepath : str or None
            Path to json file for exporting. 
            
        Returns
        -------
        None or exported Object
        """
        
        out_dict = {'populationSize': self.popSize,
                    'targetMIproportion': self.targetMIproportion,
                    'targetAUPRC':self.targetAUPRC,
                    'alpha': self.alpha,
                    'beta': self.beta,
                    }
        
        #Fills outputted orgList with an export dict for each organism
        out_dict['populationList'] = [self.populationList[i]._export() for i in range(self.popSize)]
        
        #Outputs
        with open(filepath, 'w') as f:
            json.dump(out_dict, f)
            
        
        
        
        
    

                    
            
            
if __name__ == "__main__":
    """
    Generates a sample coevolutionary population
    """
    exampleCoevPop = coevPopulation(2,0.50,1,10,2046,3,5,32,["org1","org2"])
    #exampleCoevPop = coevPopulation(1,0.50,1,10,2046,32,8,32)
    print(exampleCoevPop.runOneGeneration())
    
        
        
            
        
            
        
            
        
        
        
        
    
        