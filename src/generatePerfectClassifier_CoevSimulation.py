# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 22:40:21 2024

@author: thele
"""
from src.Objects.coev_organism import coevPopulation
import matplotlib.pyplot as plt
import json
import pickle

if __name__ == "__main__":
    #Opens the json file and reads
    with open("./settings_generatePerfectClassifier.json") as jfile:
        json_file = json.load(jfile)
        
    #Input File to run simulation on pregenerated population (Used to resume coevolutionary simulation or change parameters. None for random initialization)
    input_file = json_file["Input_File"]
    
    #Number of organisms for coevolutionary population
    num_org = json_file["Number_Organisms"]
    
    #Percent MI requested for recognizer
    percent_MI = json_file["Percent_MI"]
    
    #Parameters for fitness functions
    #Alpha - Weighting on Percent MI (1 = Equal weight as AUPRC)
    #Beta - Base of exponent for computing num offspring
    alpha = json_file["Alpha"]
    beta = json_file["Beta"]
    
    #Settings for each coevolutionary organisms. 
    genomeLength = json_file["Genome_length"]
    numBS_Recognizer = json_file["Number_BindingSites_Recognizer"]
    lenBS = json_file["Length_BindingSites"]
    numBS = json_file["Number_BindingSites"]
    
    output_file = json_file["Output_File"]
    
    #Ensures output file is specified
    if output_file is None:
        raise ValueError("Make sure to input a path for output file to output results of coevolutionary simulation")
    
    #Lists for plotting
    generationList = []
    meanAUPRCList = []
    MIList = []
    maxAUPRCList = []
    i = 0
    
    #Initializes population
    if input_file is None:
        inputPop = coevPopulation(num_org, percent_MI, alpha, beta, genomeLength, numBS_Recognizer, lenBS, numBS) 
    else:
        inputPop = coevPopulation(importFile = input_file)
    
    #Initialized for first iteration
    generationSum = {"Top AUPRC": 0}
    
    #Runs iterations of coevolutionary simulation until one organism is at 1 AUPRC
    while generationSum["Top AUPRC"] < 1:
        
        #Generates next generation and saves summary statistics
        generationSum = inputPop.runOneGeneration()
        
        generationList.append(i)
        meanAUPRCList.append(generationSum["Mean AUPRC"])
        maxAUPRCList.append(generationSum["Top AUPRC"])
        mean_IC, mean_MI, mean_TI = generationSum["Mean Information"]
        MIList.append(mean_MI/mean_TI)
        
        #Plots every 5 generations
        if i % 5 == 0:
            plt.plot(generationList, meanAUPRCList, color = 'blue', label = 'Average AUPRC')
            plt.plot(generationList,MIList, color = 'orange', label = 'Average Percent MI')
            plt.plot(generationList,maxAUPRCList, color = 'red', label = "Max AUPRC in population")
            plt.xlabel("Generations of Simulation")
            plt.legend()
            
            plt.pause(0.000001)
            
        
        print(
            f"Generation {i}: Mean Fitness - {generationSum['Mean Fitness']}, Mean % MI - {MIList[i]}, Mean AUPRC - {meanAUPRCList[i]}")


        i += 1
    
    #Grabs 1AUPRC organism
    for i in range(len(inputPop.populationList)):
        if (inputPop.populationList[i].returnAUPRC() == 1):
            topOrganism = input.populationList[i]
            

    #Outputs top organism to output file
    topOrganism._export(output_file)