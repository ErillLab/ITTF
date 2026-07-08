# -*- coding: utf-8 -*-

import pickle
import csv
import os
from scipy.stats import sem
from statistics import mean
import math
import dipwmsearch as ds
from math import log2
from scipy.spatial.distance import hamming
import numpy as np
import pandas as pd

def reachedAsymptote(yList):
    """
    Accepts 2D lists where each successive element (list) is the next generation of mutations
    
    yxList - 1D List where the first element of each list is the x value and the second element is the y value
                                                                  
    Called by loop in the main function used to induce mutation in the population. Returns true when slope hits threshold close to 0 (Asymptote).
    Returns false is slope does not hit threshold or 10 generations have not yet been generated
                                                                  
    """
    
    #Edge case
    if (len(yList) <= 30):
        return False
    
    #Computes Average of last 10 elements of list and the 10 before it
    averageYValcurr = mean([x for x in yList[-10:]])    
    averageYvalprev = mean([x for x in yList[-30:-19]])
    
    #Computes the magnitude of the difference
    differenceMagnitude = abs(averageYValcurr - averageYvalprev)
    print(differenceMagnitude)
    #Returns true if reached asymptote (slope close to 0)
    if (differenceMagnitude <= 0.0005):
        #Continues iteration until number of generations is a multiple of 100
        if ( ((len(yList) - 1)%100) == 0):
            return True
        else:
            return False
    else:
        return False
    

def equalizeGenerations(decayFile, csvFile, numGenerations):
    """
    Helper function to take a decayFile associated with the numbers in the csvFile and further decay the population till numGenerations

    Parameters
    ----------
    decayFile : string
        pickle file path containing decayed population
    csvFile : string
        csv file path containing metrics from pickle population in decayFile
    numGenerations : int
        number of generations to decay till

    Returns
    -------
    None.

    """
    
    #Loads population
    with open(decayFile,'rb') as populationFileVar:
        population = pickle.load(populationFileVar)
        
        populationFileVar.close()
    
    #Recreates lossList
    lossList = []
    
    #Loads csv file
    with open(csvFile, "r", newline = "") as csvFileHandle:
        csvreader = csv.reader(csvFileHandle)
    
        for row in csvreader:
            lossList.append(row)
            
        csvFileHandle.close()
        
    #Gets last generation number
    latestGeneration = int(lossList[-1][0])
    
    #Induces mutation until we get to same generation number as numGenerations
    while latestGeneration < numGenerations:
        
        #Mutates organism lines, updating Average IC, MI, and TI
        population.Lose_Information()
        Average_IC = sum([obj.IC for obj in population.population_list])/len(population.population_list)
        SEM_IC = sem([obj.IC for obj in population.population_list])
        
        Average_MI = sum([obj.MI for obj in population.population_list])/len(population.population_list)
        SEM_MI = sem([obj.MI for obj in population.population_list])
        
        Average_TI = sum([obj.TI for obj in population.population_list])/len(population.population_list)
        
        #Adds to generation counter
        latestGeneration += 1
        
        #Adds to list
        lossList.append([latestGeneration, Average_IC, Average_MI, Average_TI, SEM_IC, SEM_MI])
        
    
    #After decay, resaves csv file
    with open(csvFile, "w", newline = "") as outputcsv:
        csvwriter = csv.writer(outputcsv)
                
        #Writes list
        csvwriter.writerows(lossList)
        
        outputcsv.close()
    
def equalizeFiles(filePairs):
    """
    Function passed a list of file pairs (1st = decay pickle population file and 2nd = CSVfile). 
    This function determines the maxmium generation among all of the csv files and induces additional 
    mutations so that all csv files have the same # of generations

    Parameters
    ----------
    filepairs: List of tuples
        Collection of tuples where each tuple contains #1 the decayed population file path and #2 the csv file path
        

    Returns
    -------
    None.

    """    
    
    #Stores csv files in the following list
    csvFileList = []
    
    #Variable for max generation number
    maxGeneration = 0
    
    #Populates list    
    for i in range(len(filePairs)):
        csvFileList.append(filePairs[i][1])
    
    #Iterates through list looking for largest num of generations in csv file
    for filePath in csvFileList:
        with open(filePath, "r", newline = "") as csvFileHandle:
            csvreader = csv.reader(csvFileHandle)
            
            csvGeneration = int(list(csvreader)[-1][0])
            #Grabs last row and checks if the value is greater than previous maxGeneration
            if (csvGeneration) > maxGeneration:
                maxGeneration = csvGeneration
            csvFileHandle.close()
            
    #Runs equalizeGenerations on each tuple of files
    for i in range(len(filePairs)):
        equalizeGenerations(filePairs[i][0], filePairs[i][1], maxGeneration)
        
def equalizeSubFiles(DecayPopulationDir, csvFolderDir):
    """
    Runs equalizeFiles on all populations contained 

    Parameters
    ----------
    csvFolderDir : TYPE
        DESCRIPTION.
    DecayPopulationDir : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    filesList = []
    for subdir, dirs, files in os.walk(csvFolderDir):
        for file in files:
            noExtensionFile = os.path.splitext(file)[0]
            filesList.append((DecayPopulationDir + noExtensionFile + ".pickle", csvFolderDir + noExtensionFile + ".csv"))
    
    equalizeFiles(filesList)
            
def CalculateJennsenShannon(TrueBS, DivergedBS, N):
    '''
    Function to determine the divergence between two probability distribution (two lists
    of binding sites: TrueBS and DivergedBS). 
    
    This function computes the frequency of each N-mer in each list of binding sites.
    The divergence is then computed through a weighted log-likelihood ratio (JS Divergence).

    Note: this function is restricted to 1-mers or 2-mers to mitigate information overestimation. 
    
        
    Parameters
    ----------
    TrueBS : List
        List of TFBS of True Motif
    DivergedBS : List
        List of TFBS of Diverged Motif
    N : Int
        N-mer comparison (N = 1 or 2). Used to grab Nmers and calculate frequency

    Returns
    -------
    JSD

    '''
    #Initializes true and diverged frequency dictionary
    trueFreqDict = {}
    divergedFreqDict = {}
    mixedDict = {}
    
    alphabet = ["A","C","G","T"]
    alphabet2 = ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT",
                 "TA","TC","TG","TT"]
    
    #Limits runs to motifs of length 1 or 2
    if ((len(TrueBS[0]) <= 2) & (len(DivergedBS[0]) <= 2)):
        print(TrueBS[0], DivergedBS[0])
        if N == 1:
            for letter in alphabet:
                trueFreqDict[letter] = 0
                divergedFreqDict[letter] = 0
                mixedDict[letter] = 0
    
            
        elif N == 2:
            for letter1 in alphabet:
                for letter2 in alphabet:
                    trueFreqDict[letter1 + letter2] = 0
                    divergedFreqDict[letter1 + letter2] = 0
                    mixedDict[letter1 + letter2] = 0
        
        
        
        #Populates trueFreqDict with counts
        for bindingSite in TrueBS:
            #Adds N-mer to count
            if N == 1:
                trueFreqDict[bindingSite[0]] += 1
            elif N==2:
                trueFreqDict[bindingSite[0] + bindingSite[1]] += 1
        
        #Populates divergedFreqDict with counts
        for bindingSiteD in DivergedBS:
            #Adds N-mer to count
            if N == 1:
                divergedFreqDict[bindingSiteD[0]] += 1
            elif N==2:
                divergedFreqDict[bindingSiteD[0] + bindingSiteD[1]] += 1
                
        #Converts to frequency 
        sumTrueCounts = sum(trueFreqDict.values())
        sumDivergedCounts = sum(divergedFreqDict.values())
        
        #Changes true counts to frequencies
        for i in trueFreqDict.keys():
            trueFreqDict[i] = trueFreqDict[i]/sumTrueCounts
    
        #Changes diverged counts to frequencies
        for i in divergedFreqDict.keys():
            divergedFreqDict[i] = divergedFreqDict[i]/sumDivergedCounts
        
        #Generates mixed distribution
        for i in mixedDict.keys():
            mixedDict[i] = (trueFreqDict[i] + divergedFreqDict[i])/2
       
        #Divergence formula used below (weighted log likelihood ratio)
        Divergence = lambda P, Q: P* log2(P/Q) 
        
        if N==1:
            #Initializes KLD
            KLD_P = 0
            KLD_Q = 0
            for state in alphabet:
                               
                #Replaces 0 with integer value if P(x) is not zero for state (Otherwise adds 0)
                if (trueFreqDict[state] != 0):
                    KLD_P += Divergence(trueFreqDict[state], mixedDict[state])
                    
                #Replaces 0 with integer value if Q(x) is not zero for state (Otherwise adds 0)
                if (divergedFreqDict[state] != 0):
                    KLD_Q += Divergence(divergedFreqDict[state], mixedDict[state])
                
            
        
        if N ==2:
            
            #Initializes KLD
            KLD_P = 0
            KLD_Q = 0
            
            #Computes KLD
            for state in alphabet2:
                                
                #Replaces 0 with integer value if P(x) is not zero for state (Otherwise adds 0)
                if (trueFreqDict[state] != 0):
                    KLD_P += Divergence(trueFreqDict[state], mixedDict[state])
                    
                #Replaces 0 with integer value if Q(x) is not zero for state (Otherwise adds 0)
                if (divergedFreqDict[state] != 0):
                    KLD_Q += Divergence(divergedFreqDict[state], mixedDict[state])
                
        
        #Returns JSD as the average KLD
        return ((KLD_P + KLD_Q)/2)
    
    else:
        return 1;

def AverageJSD(initialTFBS, newTFBS, N):
    """
    Function used in population to Compute JSD between correct set of motifs and return list describing divergence of each organism
    in the population from initialTFBS

    Parameters
    ----------
    initialTFBS : TYPE
        DESCRIPTION.
    newTFBS : TYPE
        DESCRIPTION.
    N : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    JSDList = []
    for i in range(len(newTFBS)):
        JSDList.append(CalculateJennsenShannon(initialTFBS[i], newTFBS[i], N))
        print(CalculateJennsenShannon(initialTFBS[i], newTFBS[i], N))
    return JSDList

def create_dipwm_list(TFBS_list):
    """
    Helper Function that accepts Transcription Factor Binding sites, computes dinucleotide
    frequencies at each position (16 potential states), and creates diNucleotide
    list appropriate for creating Kulaokovsky diPWM object
    
    This function does this by populating a list of dictionaries with the dinucleotide counts 
    at each position of the motif, converting counts to scores using the diChipMunk PWM formula,
    and writing to a .dpwm file path provided to the function in the parameters (The format for these files 
    entail a position x dinucleotide matrix separated by tabs)
    
    diChIPMunk formula: log(x + cq/((N + c)q)) where x = counts of a dinucleotide at a position, c = pseudocount (ln(N)),
    q = background frequency of nucleotide (uniform), N = number of sequences
    Parameters
    ----------
    TFBS_list: List
        List of Transcription Factor Binding sites.
    
        
    Returns
    -------
    scoreList: list
        List of dictionaries with dinucleotide frequencies at each position. 
        Used by coevolutionary Objects

    """
    
    #Initializes alphabet (16 dinucleotides)
    diNucAlphabet = ["AA", "AC", "AG", "AT", "CA", "CC", 
                     "CG", "CT", "GA", "GC", "GG", "GT",
                     "TA", "TC", "TG", "TT"]
    
    #Initailizes frequency list which contains dictionaries of diNuc freqs for each position
    diNucFreqs = []
    
    #Initializes variable holding length of binding sites
    lenBS = len(TFBS_list[0])
    
    #Populates with Counts
    for i in range(lenBS-1):
        newDict = {}
        for diNuc in diNucAlphabet:
            newDict[diNuc] = 0
        diNucFreqs.append(newDict)
    
    #Fills diNucFreqs with counts at each position
    for i in range(0, lenBS-1):
        
        #Iterates through each binding site in list
        for bsNum in range(0,len(TFBS_list)):
            site = TFBS_list[bsNum]
            
            #Adds count to diNucFreq
            diNucFreqs[i][site[i] + site[i+1]] = diNucFreqs[i][site[i] + site[i + 1]] + 1
    
    #Variables related to computing score in diChipMunk manner
    c_Munk = math.log(len(TFBS_list)) #Pseudocount parameter
    q_Munk = 1/16 #Background frequency (uniform)
    N_Munk = len(TFBS_list) #Number of sequences
    
    
    #Converts to scores
    for i in range(len(diNucFreqs)):
        
        #Iterates across alphabet to compute score at each position
        for diNuc in diNucAlphabet:
            diNucFreqs[i][diNuc] = math.log(
                (diNucFreqs[i][diNuc] + (c_Munk * q_Munk))/((N_Munk + c_Munk) * q_Munk)
                                            ,2)
            
            
    
    #Begins to write to file by creating lists for the scores of each dinuc at each positon
    scoreList = []
    
    #Iterates across position
    for i in range(len(diNucFreqs)):
        #Adds potential scores at a position to a list and puts it in scoreList
        positionScore = [score for score in diNucFreqs[i].values()]
        scoreList.append(positionScore)
    
    #returns file
    return scoreList

def create_dipwm_object(TFBS_list):
    """
    Function used to convert list of transcription factor binding sites
    into a diPWM object as defined by the Kulakovsky package
    
    Used in coevolutionary algorithm and genome search experiment
    to compute scores of recognizer on each position of the genome
    
    Parameters
    ----------
    scoreList : TYPE
        DESCRIPTION.

    Returns
    -------
    Kulakovsky diPWM object

    """
    scoreList = create_dipwm_list(TFBS_list)
    return ds.diPWM(scoreList)

def averageHamming(coevPop, geneLength):
    geneList = []
    HammingList = []
    for i in coevPop.populationList:
        geneList.append(i.genomeObj.genome[:geneLength])
    for i in range(len(geneList)-1):
        tempList = []
        
        for j in range(i+1,len(geneList)):
            tempList.append(hamming(list(geneList[i]), list(geneList[j])))
        HammingList.append(mean(tempList))
    
    return mean(HammingList)


def compute_seq_logo_matrix(seqs):
    '''
    Returns the matrix containing the data for the "sequence logo" generated from
    `seqs`, which is the list of (aligned) sequences. The matrix is returned as
    a (pandas) dataframe (Lx4, where L is the length of each sequence), which
    contains the height of each letter at each position of the sequence logo.
    
    Sequences are assumed to be strings of DNA, containing only the letters
    'A', 'C', 'G' and 'T' (uppercase), all the same length.
    '''
    N = len(seqs)
    alphabet = ['A', 'C', 'G', 'T']
    mat = []
    for pos in range(len(seqs[0])):
        counts = {nt: 0 for nt in alphabet}
        for s in seqs:
            counts[s[pos]] += 1
        freqs = {nt: counts[nt] / N for nt in alphabet}
        # Shannon entropy
        H = -sum(p * np.log2(p) for p in freqs.values() if p > 0)
        SI = 2 - H
        # Height of each letter at this position of the sequence logo
        mat.append({nt: freqs[nt] * SI for nt in alphabet})
    return pd.DataFrame(mat)



