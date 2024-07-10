# -*- coding: utf-8 
import json
import numpy as np
from math import log2
import csv
import pickle
from Objects.population import Population
from Bio import SeqIO
from scipy.stats import sem
from statistics import mean 
from copy import deepcopy
from helperFunctions import reachedAsymptote

def ComputeKL(Px, Qx):
    '''
    Computes KL as P || Q
    Returns value

    Parameters
    ----------
    Px : TYPE
        DESCRIPTION.
    Qx : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    return Px * log2((Px)/(Qx))

def CalculateJennsenShannon(TrueBS, DivergedBS, N):
    '''
    Function to determine the divergence between the probability distributions of 
    TrueBS and DivergedBS given that TrueBS is the true distribution. Here, 
    we compare the frequency of N-mers, and we restrict our analysis to 
    1-mers or 2-mers to mitigate information overestimation. This function 
    permits us to quantify the divergence of a query binding site to the 
    TF's binding sites, which describe its binding tendency.
    
    KLD is computed by generating a log-likelihood for each x in two probability 
    distribution and taking the sum of the log-likelihood ratios weighted by the 
    frequency of x in the true probability distribution. As such, KLD represents
    the average log-likelihood ratio of two distributions given the first one 
    is the true distribution. 
    
    We do not simply use KLD, however. Due to problems with zero binning, we use the Jennsen Shannon Divergence:
    The average KLD between each distribution and its mixture distribution (the distributed computed as the average of probability density of both distributuions for each state)
        
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
    KLD

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
       

        if N==1:
            #Initializes KLD
            KLD_P = 0
            KLD_Q = 0
            for state in alphabet:
                               
                #Replaces 0 with integer value if P(x) is not zero for state (Otherwise adds 0)
                if (trueFreqDict[state] != 0):
                    KLD_P += ComputeKL(trueFreqDict[state], mixedDict[state])
                    
                #Replaces 0 with integer value if Q(x) is not zero for state (Otherwise adds 0)
                if (divergedFreqDict[state] != 0):
                    KLD_Q += ComputeKL(divergedFreqDict[state], mixedDict[state])
                
            
        
        if N ==2:
            
            #Initializes KLD
            KLD_P = 0
            KLD_Q = 0
            
            #Computes KLD
            for state in alphabet2:
                                
                #Replaces 0 with integer value if P(x) is not zero for state (Otherwise adds 0)
                if (trueFreqDict[state] != 0):
                    KLD_P += ComputeKL(trueFreqDict[state], mixedDict[state])
                    
                #Replaces 0 with integer value if Q(x) is not zero for state (Otherwise adds 0)
                if (divergedFreqDict[state] != 0):
                    KLD_Q += ComputeKL(divergedFreqDict[state], mixedDict[state])
                
        
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
if __name__ == "__main__":
    '''
    Json parameters:
        populationFile (Pickle object containing path to pregenerated population)
        num_organisms - Int 
        target_TI - Int representing the total TI the organisms should attempt to reach for
        Prop_IC - Prop representing how much of this TI should be for IC, allowing the fitness function to select for organisms of specific IC/MI strategies
        MSE_thres - Float representing average error of IC and MI permitted for an organism line to have met the requirement
        populationOutputFile - Path to file to which the population will be dumped at
        propTI_terminate - Threshold of proportion of TI lost before stopping mutation
    '''
    
    #Opens the json file and reads
    with open("./settings_divergence.json") as jfile:
        json_file = json.load(jfile)
        
    #Loads file with population object 
    populationFile = json_file["Population_File"]
    
    #Loads json input into class variables
    num_organisms = json_file["Number_Organisms"]
    num_BS = json_file["Number_BindingSites"]
    length_BS = json_file["Length_BindingSites"]
    gc_cont = json_file["GC_cont"]
    mut_rate = json_file["Mutation_Rate"]
    read_in = json_file["ReadInFiles"]
    MI_mode = json_file["MI_mode"]
    
    #Loads in parameters relavent to the fitness function
    target_TI_Column = json_file["Target_TI_Column"]
    target_TI = target_TI_Column * length_BS
    Goal_IC = json_file["Prop_IC"] * target_TI
    Goal_MI = (1 - json_file["Prop_IC"]) * target_TI

    #Loads json input dealing with MSE threshold and output File creation
    MSE_thres = json_file["Eucledian_Thres"]
    populationOutputFile = json_file["populationOutputFile"]
    
    #Loads json input dealing with loss of information of population and num_generations in output file
    numGenerationsPrint = json_file["NumGenerationPrint"]
    numGenerationsTerminate = json_file["NumGenerationTerminate"]
    outputCSV = json_file["lossOutputCSV"]
    
    
    if populationFile == "":
        
        #Initailizes population object with parameters
        population = Population(num_organisms,num_BS,length_BS,gc_cont,mut_rate, Goal_IC, Goal_MI, constrictMode = "", MI_mode_org = MI_mode, read_in_org_list = read_in)
        
        aboveThreshold = True
        
        #Repeats reproduce and select until metThreshold is true (when MSE values for each organism line is below threshold)
        while (aboveThreshold):
            
            #If provided population File, skips generation of population
    
                
            #Sets threshold value to false (Will iterate through MSE values to find an above threshold case)
            aboveThreshold = False
            
            #Reproduces and selects 
            population.Reproduce()           
            population.Select()
            
            
            for i in range(len(population.population_list)):
                if (population.MSE_list[i] > MSE_thres):
                    aboveThreshold = True
                    print(round(population.MSE_list[i],4))

        
        if populationOutputFile != "":
            with open(populationOutputFile, 'wb') as output:
                pickle.dump(population, output)
                output.close()
            
    else:
        #Else loads population file (after selection)
        with open(populationFile,'rb') as populationFileVar:
            population = pickle.load(populationFileVar)
    
    #Prepares list containing sublists with the elements in fields
    lossList = []
    abrevLossList = []
    frozenMotif = deepcopy([obj.TFBS for obj in population.population_list])
    fields = ["# Generation after mutation", "JSD", "SEM_JSD"]
        
    #Counter variable for appending to lossList
    generationNum = 0
    
    #Runs until average loss of TI meets threshold
    #while generationNum < numGenerationsTerminate:
    while (not reachedAsymptote(abrevLossList)):   
        #Depending on input in json file, adds to lossList
        if ((generationNum % numGenerationsPrint) == 0):
            newTFBS = [obj.TFBS for obj in population.population_list]
            JSDList = AverageJSD(frozenMotif,newTFBS,2)
            lossList.append([generationNum, mean(JSDList), sem(JSDList)])
            abrevLossList.append(mean(JSDList))
            
        #Mutates organism lines, updating Average IC, MI, and TI
        population.Lose_Information2()
        
        #Adds to generationNum
        generationNum += 1
        
    
    #Prints to csv file
    with open(outputCSV, "w") as outputcsv:
        csvwriter = csv.writer(outputcsv)
        
        #Write fields
        csvwriter.writerow(fields)
        
        #Writes list
        csvwriter.writerows(lossList)
        

    
            
    
    
        
        
    #Generates output variable if populationFile is not given
            
        
        