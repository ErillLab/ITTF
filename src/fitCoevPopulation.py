# -*- coding: utf-8 -*-
"""
Script used to fit an inputted coevPopulation (random initialization, passed in by file, etc.)
by running the simulation and selecting for offspring with recognizers that better classify their binding sites

- Used to generate a 1AUPRC organism to run the competition experiment (runCompetition.py)
- Used to run drift experiment by inputting a file path to a coevOrganism (or coevPopulation)
and the # of generations to run the experiment. 
"""
from src.Objects.coev_organism import coevOrganism, coevPopulation
import matplotlib.pyplot as plt
import json
import csv

def finishedIteration(genNum, numTerminate, topAUPRC):
    """
    Helper function to determine when to end the simulation. If numTerminate is specified,
    returns True if genNum > numTerminate. Otherwise, returns True if topAUPRC ==1
    """
    if numTerminate:
        return True if genNum > numTerminate else False
    else:
        return True if topAUPRC >= 1 else False
    
if __name__ == "__main__":
    """
    JSON_File parameters:
        Input_Files: list, str, or None
            Set to list to initialize population from list of coev_organism paths. Set to str 
            to initialize population from file path. Set to none to randomly initialize 
            population
        
        Number_Organisms: int
            Number of organisms in coevolutionary simulation population. Duplicates
            organisms in Input_files.
        
        ... (Parameters for coevPopulation. For a more detailed overview, take a look at src/Objects/coev_organisms)
        
        Beta: float
            Base of exponent for # of offspring. Affects selection vs drift of algorithm
        
        Number_Generations: int or None
            Number of generations to run simulation for. 
            If none, runs until 1AUPRC organism is obtained. 
        
        Output_CSV_File: str or None
            Path to save csv file to with statisticsc about the population after
            each generation of simulation (Max AUPRC, Mean AUPRC, Mean MI). 
            Set to none to not save 
            
        Ouput_Org_File: str or None
            Path to save 1AUPRC organism to (JSON). If 1AUPRC organism has not been found, 
            outputs top organism. Set Num_Generations to none to guarantee 1AUPRC organism is found.
            Set to none to not save top organism
        
        Output_Pop_File: str or None
            Path for output population (JSON). Set to none to not save
            
    """
    
    #Opens the json file and reads
    with open("./src/settings_fitCoevPop.json") as jfile:
        json_file = json.load(jfile)
        
    #Input File to run simulation on pregenerated population (Used to resume coevolutionary simulation or change parameters. None for random initialization)
    input_files = json_file["Input_Files"]
    
    #Number of organisms for coevolutionary population
    num_org = json_file["Number_Organisms"]
    
    #Percent MI requested for recognizer
    percent_MI = json_file["Percent_MI"]
    
    #Parameters for fitness functions
    #Alpha - Weighting on Percent MI (1 = Equal weight as AUPRC)
    #Beta - Base of exponent for computing num offspring
    alpha = json_file["Alpha"]
    beta = json_file["Beta"]
    
    #Settings for each coevolutionary organism. Must only be specified if input_files is None
    genomeLength = json_file["Genome_Length"]
    numBS_Recognizer = json_file["Number_BindingSites_Recognizer"]
    lenBS = json_file["Length_BindingSites"]
    numBS = json_file["Number_BindingSites"]
    
    num_gen = json_file["Number_Generations"]
    
    #Grabs output files  
    output_csv_file = json_file["Output_CSV_File"]
    output_org_file = json_file["Output_Org_File"]
    output_pop_file = json_file["Output_Pop_File"]
    
    #Lists for plotting
    generationList = []
    meanAUPRCList = []
    averageMIList = []
    maxAUPRCList = []
    i = 0
    
    #Initializes population from list of organisms
    if type(input_files) is list:
        inputPop = coevPopulation(num_org, percent_MI, alpha, beta, genomeLength, numBS_Recognizer, lenBS, numBS,importList = [coevOrganism(importJSON = filePath) for filePath in input_files])
    
    #Initializes population from file path
    elif type(input_files) is str:
        inputPop = coevPopulation(importFile = input_files)
    
    #Random Initialization
    else:
        inputPop = coevPopulation(num_org, percent_MI, alpha, beta, genomeLength, numBS_Recognizer, lenBS, numBS) 
    
    #Initialized for first iteration
    generationSum = {"Top AUPRC": 0}
    
    #Runs iterations of coevolutionary simulation
    while not finishedIteration(i,num_gen,generationSum["Top AUPRC"]):
        
        #Generates next generation and saves summary statistics
        generationSum = inputPop.runOneGeneration()
        
        generationList.append(i+1)
        meanAUPRCList.append(generationSum["Mean AUPRC"])
        maxAUPRCList.append(generationSum["Top AUPRC"])
        mean_IC, mean_MI, mean_TI = generationSum["Mean Information"]
        averageMIList.append(mean_MI/mean_TI)
        
        #Plots every 5 generations
        if i % 5 == 0:
            plt.plot(generationList, meanAUPRCList, color = 'blue', label = 'Average AUPRC')
            plt.plot(generationList,averageMIList, color = 'orange', label = 'Average Percent MI')
            plt.plot(generationList,maxAUPRCList, color = 'red', label = "Max AUPRC in population")
            plt.xlabel("Generations of Simulation")
            plt.legend()
            
            plt.pause(0.000001)
            
        
        print(
            f"Generation {i+1}: Mean Fitness - {generationSum['Mean Fitness']}, Mean % MI - {averageMIList[i]}, Mean AUPRC - {meanAUPRCList[i]}")


        i += 1
    
    #Saves CSV file
    if output_csv_file:
        with open(output_csv_file, mode = "w", newline = '') as output:
            csvwriter = csv.writer(output)
            csvwriter.writerow(["Generation", "Mean AUPRC", "Max AUPRC", "Average MI"]) #Writes names of columns
            
            #Writes info to csv file 
            for i in range(len(generationList)):
                csvwriter.writerow([generationList[i],meanAUPRCList[i], maxAUPRCList[i], averageMIList[i]])
            output.close()
        
    #Grabs Best classifier and saves to file 
    for i in range(len(inputPop.populationList)):
        if (inputPop.populationList[i].returnAUPRC() == generationSum["Top AUPRC"]):
            topOrganism = inputPop.populationList[i]
    
    if output_org_file:
        topOrganism._export(output_org_file)
    
    #Saves population file
    if output_pop_file:
        inputPop._export(output_pop_file)