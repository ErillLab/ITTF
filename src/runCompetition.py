"""
Runs competition experiment with 1AUPRC organisms generated from fitCoevPopulation.py.
This entails loading in MI/IC organisms into a population, running the simulation
for a set number of generations, and outputting the frequency of MI/IC organisms
to a csv file

"""
import pickle
import csv
import copy
import matplotlib.pyplot as plt
import json
from src.Objects.coev_organism import coevPopulation, coevOrganism

if __name__ == "__main__":
    """
    JSON_File parameters:
        Input_Files: list
            List of file paths for each 1AUPC coev_organism. JSON file paths
            
        ID_list: list
            list of ID strings for each ancestral organism in Input_Files. IDs are inherited by offspring
            
        Number_Organisms: int
            Number of organisms in coevolutionary simulation population. 
            if Number_Organisms > len(Input_Files), then organisms are duplicated to fill population
            
        Beta: float
            Base of exponent for # of offspring. Affects selection vs drift of algorithm
        
        Number_Generations: int
            Number of generations to run simulation for
            
        Ouput_File: str
            Path for output file (csv). Contains statistics about population and 
            frequency of MI/IC organisms
        
    """
    with open("./src/settings_runCompetition.json", 'rb') as jfile:
        json_file = json.load(jfile)
    
    #Input File to run simulation on pregenerated population 
    input_files = json_file["Input_Files"]
    
    #ID_list. Must be as long as input_files
    id_list = json_file["ID_list"]
    if (len(id_list) != len(input_files)):
        raise ValueError("# of IDs does not equal # of input coevolutionary organisms")
    
    #Number of organisms for coevolutionary population
    num_org = json_file["Number_Organisms"]
    
    #Number of generations to run simulation
    num_gen = json_file["Number_Generations"]
    
    #Parameter for fitness functions
    #Beta - Base of exponent for computing num offspring
    beta = json_file["Beta"]
    
    #Output file
    output_file = json_file["Output_File"]
    
    #Grabs orgs
    input_orgs = []
    
    #Adds organisms to list with IDs
    for i in range(len(input_files)):
        importOrg = coevOrganism(importJSON = input_files[i])
        importOrg.orgID = id_list[i]
        input_orgs.append(importOrg)
            
    
    motifPop = coevPopulation(popSize = num_org, targetMIproportion = 0.5, 
                              alpha = 0,beta=beta, importList = input_orgs)

    #Keeps track of results from each generation of coevolutionary simulation
    generationList = []
    meanAUPRCList = []
    maxAUPRCList = []
    averageMIList = []
    freqOrgsList = {tag:[] for tag in id_list}
    
    
    #Runs simulation
    for n in range(num_gen):
        #Runs one generation
        generationSum = motifPop.runOneGeneration()
        
        #Adds to lists for mean AUPRC, max AUPRC, and mean information
        generationList.append(n)
        meanAUPRCList.append(generationSum["Mean AUPRC"])
        maxAUPRCList.append(generationSum["Top AUPRC"])
        averageMIList.append(generationSum["Mean Information"][1]/generationSum["Mean Information"][2])
                   
        #Prepares dictionary that will contain counts of organisms with each ID
        freqOrgsID = {tag:0 for tag in id_list}
        
        #Iterates through populationList
        for i in range(motifPop.popSize):
            referenceOrg = motifPop.populationList[i]
            
            #Adds count to appropriate key
            for tag in freqOrgsID.keys():
                if referenceOrg.orgID == tag:
                    freqOrgsID[tag] += 1
                    
        #Converts to frequency
        freqOrgsID = {tag:freqOrgsID[tag]/num_org for tag in id_list}
                   
        #Adds to list
        for tag in id_list:
            freqOrgsList[tag].append(freqOrgsID[tag])
                    
        #avgHamList.append(averageHamming(testPop, 4*12))
        if n % 5 == 0:
            #plt.plot(generationList, AUPRCList, color = 'blue', label = 'AUPRC')
            plt.plot(generationList,averageMIList, color = 'orange', label = 'Average Percent MI of entire population')
            
            #Plots frequency of first tag
            plt.plot(generationList,freqOrgsList[id_list[0]], color = "purple", label = "Percent of population from " + id_list[0] + " organism")

            plt.xlabel("Generations of Simulation")
            plt.legend(loc = "upper left", bbox_to_anchor=(1,1))
            plt.ylim(0,1)
            plt.pause(0.000001)
        print(f"Generation {n}") 
    
    #Writes AUClist to file
    with open(output_file, mode = "w", newline = '') as output:
        csvwriter = csv.writer(output)
        csvwriter.writerow(["Generation", "Mean AUPRC", "Max AUPRC", "Average MI"] + [f"{tag} org frequency" for tag in id_list]) #Writes names of columns
        
        #Writes info to csv file 
        for i in range(len(generationList)):
            csvwriter.writerow([generationList[i],meanAUPRCList[i], maxAUPRCList[i], averageMIList[i]] + 
                               [freqOrgsList[tag][i] for tag in id_list])
        output.close()
