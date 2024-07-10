# -*- coding: utf-8 -*-
"""
Genome Search Experiment
@author: thele
"""

import csv 
import math
from math import log2
import json
import dipwmsearch as ds
import pickle
from sklearn import metrics
from Objects.population import Population
from Objects.genome import Genome, CreateDiPWMfile, scoreGenome
        

# XXX

if __name__ == "__main__":
    #Opens the json file and reads
    with open("./settings_GenomeExperiment.json") as jfile:
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
    
    #Possible values: "", "IC", "MI"
    constrictMode = json_file["ConstrictMode"]
    if ((constrictMode != "IC") & (constrictMode != "MI")):
        constrictMode = ""
        
        

    #Loads json input dealing with MSE threshold and output File creation
    MSE_thres = json_file["Eucledian_Thres"]
    populationOutputFile = json_file["populationOutputFile"]
    
    #Loads json input dealing with loss of information of population and num_generations in output file
    numGenerationsPrint = json_file["NumGenerationPrint"]
    numGenerationsTerminate = json_file["NumGenerationTerminate"]
    outputCSV = json_file["lossOutputCSV"]
    
    #Computes Genome Length based on Schenider's equation R_freq = R_seq
    R_seq = target_TI
    GenomeLength = 2**(R_seq + log2(num_BS)) #R_seq = R_freq = log2(G) - log2(y)
    GenomeLength = math.ceil(GenomeLength)
    
    if (GenomeLength < (num_BS * length_BS)):
        print("PROBLEM!")
    
    #Genome length minus Genome positions taken by binding sites
    nonBSGenomeLength = GenomeLength - (length_BS * num_BS)
    BSGenomeLength = length_BS * num_BS
    
    if nonBSGenomeLength < 0:
        print("PROBLEM!")
    
    if populationFile == "":
        
        #Initailizes population object with parameters
        population = Population(num_organisms,num_BS,length_BS,gc_cont,mut_rate, Goal_IC, Goal_MI, constrictMode, MI_mode_org = MI_mode, read_in_org_list = read_in)
        
        aboveThreshold = True
        
        #Repeats reproduce and select until metThreshold is true (when MSE values for each organism line is below threshold)
        while (aboveThreshold):
            
            #If provided population File, skips generation of population
            
            numAboveThreshold = 0 #Variable to keep track of how many organisms are above the threshold
                
            #Sets threshold value to false (Will iterate through MSE values to find an above threshold case)
            aboveThreshold = False
            
            #Reproduces and selects 
            population.Reproduce()           
            population.Select()
            
            
            for i in range(len(population.population_list)):
                if (population.MSE_list[i] > MSE_thres):
                    aboveThreshold = True
                    numAboveThreshold += 1
             
                       
            print(numAboveThreshold, round(population.MSE_list[-1],2))

        
        if populationOutputFile != "":
            with open(populationOutputFile, 'wb') as output:
                pickle.dump(population, output)
                output.close()
            
    else:
        #Else loads population file (after selection)
        with open(populationFile,'rb') as populationFileVar:
            population = pickle.load(populationFileVar)
    
    #raise ValueError("hm")
    #Initializes AUClist, a 2D list which contains a list of AUCs for each organism at each generation of mutation  
    AUClist = []
    organismAUCListCollection = []
    NumberOrganisms = 0 #Keeps track of number of organisms
    
    organismCounter = 0
    #Assembles genome, computes PR curve, and induces mutations for each organism in population
    for organism in population.population_list:
        print(organismCounter); organismCounter +=1
        #Initializes list containing letter at each position
        genomeList = []
        bindingSites = organism.TFBS
    
        #Creates genome object from random letters and binding sites
        genome = Genome(nonBSGenomeLength, bindingSites)
        
        #Creates positive negative list describing each index of the genome as a positive or negative
        PNlist = genome.positiveNegativeList #Length of circular genome. Contains 1 for a position that is a binding site and 0 for a non binding site
        
        #Creates motif        
        CreateDiPWMfile(organism.TFBS, "", os.getcwd() + "/test.dpwm")
        
        #Reads motif
        diP = ds.create_diPwm(os.getcwd() + "/test.dpwm")
        
        #Scores all positions in genome 
        scoresList =  scoreGenome(diP, length_BS, genome.getCircularGenome())
        
        if len(scoresList) != (GenomeLength):
            raise ValueError("Problem!") 
        
                
        #The TF can be described as a classifier of positives and negatives of potential binding sites
        #We use the skilearn Precision, recall functions to quantify how good the classifer is 
        organismAUCList = []
        precisionArray, recallArray, thresholds = metrics.precision_recall_curve(PNlist, scoresList)
        
        #Computes AUC as an indicator of how good of a job the TF does
        organismAUCList.append(metrics.auc(recallArray, precisionArray))
        
        #Induces mutations on genome, recomputes scores, and recomputes PR
        for i in range(numGenerationsTerminate-1):
            
            #Induces mutation
            genome.induceMutation(mut_rate)
            
            #Computes Scores
            scoresList = scoreGenome(diP, length_BS, genome.getCircularGenome())                   
            if len(scoresList) != (GenomeLength):
                raise ValueError("Problem!")
                
            
            #The TF can be described as a classifier of positives and negatives of potential binding sites
            #We use the skilearn Precision, recall functions to quantify how good the classifer is 
            precisionArray, recallArray, thresholds = metrics.precision_recall_curve(PNlist, scoresList)
            
            #Computes AUC as an indicator of how good of a job the TF does
            organismAUCList.append(metrics.auc(recallArray, precisionArray))
            
        #Adds to 2D list
        organismAUCListCollection.append(organismAUCList)
        NumberOrganisms += 1
        
    
    #After all organisms have been genome searched, averages AUC at each generation of mutation and creates final, averaged AUC list
    for i in range(numGenerationsTerminate):
        sumAUC = 0
        
        #Iterates through each list to get AUC at specific generation
        for auclist in organismAUCListCollection:
            sumAUC += auclist[i]
        
        #Goes from sum to average
        meanAUC = sumAUC/NumberOrganisms
        
        #Appends to list
        AUClist.append(meanAUC)
        

    
    #Writes AUClist to file
    with open(outputCSV, mode = "w", newline = '') as output:
        csvwriter = csv.writer(output)
        csvwriter.writerow(["Generation", "AUC"]) #Writes names of columns
        
        #Writes generation and AUC 
        for i in range(len(AUClist)):
            csvwriter.writerow([i, AUClist[i]])
        
        
    """
    #Initializes list to be iterated through
    inc = (maximumScore-minimumScore)/199
    thresholdList = np.arange(minimumScore, maximumScore, inc)
    
    #Adds last score if rounding error did not already add it
    if len(thresholdList == 199):
    np.append(thresholdList, maximumScore)
    
    #As such, for various threshold values, we determine TP, FP, TN, and FN to compute PR (Adds pseudocount to maximumScore so that threshold includes maximumScore)
    for threshold in thresholdList:
    #Initializes variables
    TP = 0; FP = 0; TN = 0; FN = 0;
    #Iterates through scores, adding to appropriate count
    for start_position, word, score in scoresList:
        #Classified as a positive by the TF
        if round(score,7) >= round(threshold,7):
            
            #True Positive
            if start_position in trueBSIndex:
                TP += 1
            #Otherwise False Positive
            else:
                FP += 1
        #Otherwise counted as a negative by the TF
        else: 
            
            #False Negative
            if start_position in trueBSIndex:
                FN += 1
            
            #True Negative
            else:
                TN += 1
    
    #Computes Precision and Recall
    Precision =  TP/(TP+FP)
    Recall = TP/(TP+FN)
    PR_list.append((Precision,Recall))
                
    """
            

        
        
