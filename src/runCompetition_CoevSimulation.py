"""
Runs competition experiment with 1AUPRC organisms generated from generatePerfectClassifer_CoevSimulation.py
JSON_File parameters:
    -Input_Files = List of file paths for each 1AUPC coev_organism. JSON files
    -ID_list = IDs for each ancestral organism in Input_Files. IDs are maintained by offspring
    -Number_Organisms = Number of organisms in coevolutionary simulation population. Duplicates
    organisms in Input_files.
    -Beta = Base of exponent for # of offspring. Affects selection vs drift of algorithm
    -Ouput_File = Path for output file (csv)
    
"""
import pickle
import csv
import copy
import matplotlib.pyplot as plt
import json
from src.Objects.coev_organism import coevPopulation, coevOrganism

if __name__ == "__main__":
    with open("./src/settings_runCompetition.json", 'rb') as jfile:
        json_file = json.load(jfile)
    
    #Input File to run simulation on pregenerated population (Used to resume coevolutionary simulation or change parameters. Unnecessary)
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
    for i in range(input_files):
        importOrg = coevOrganism(importJSON = input_files[i])
        importOrg.orgID = id_list[i]
        input_orgs.append(importOrg)
            
    #Generates population list
    population_list = []
    id_list = []
    
    for i in range(num_org):
        #Adds equal number of each organism to populationList 
        population_list.append(copy.deepcopy(input_orgs[i%len(input_files)]))
        
    #Initializes population
    inputPop = coevPopulation(num_org,0,0,beta,input_orgs[0].genomeSize,
                              input_orgs[0].motif_numberBS,input_orgs[0].lengthBS, 
                              input_orgs[0].numberBS) 

    #Keeps track of results from each generation of coevolutionary simulation
    generationList = []
    meanAUPRCList = []
    maxAUPRCList = []
    
    averageMIList = []
    idOrgFrequencyList = [] #2D list with length = # of tags. 2D list of frequencies of each tag
    
    #Runs simulation
    for n in range(10000):
        #Runs one generation
        generationSum = inputPop.runOneGeneration()
        
        #Adds to lists for mean AUPRC, max AUPRC, and mean information
        generationList.append(n)
        meanAUPRCList.append(generationSum["Mean AUPRC"])
        maxAUPRCList.append(generationSum["Top AUPRC"])
        averageMIList.append(generationSum["Mean Information"][1]/generationSum["Mean Information"][2])
                   
        #Prepares dictionary that will contain counts of organisms with each ID
        numOrgsID = [{tag:0} for tag in id_list]
        
        #Iterates through populationList
        for i in range(inputPop.popSize):
            referenceOrg = inputPop.populationList[i]
            
            #Adds count to appropriate key
            for tag in numOrgsID.keys():
                if referenceOrg.orgID == tag:
                    numOrgsID[tag] += 1
        #Converts to frequency
        for i in range(len(numOrgsID.keys())):
            currentTag = numOrgsID.keys[i]
            numOrgsID[currentTag] = numOrgsID[currentTag]/num_org
            idOrgFrequencyList[i].append(numOrgsID[currentTag])
                   
                    
        #avgHamList.append(averageHamming(testPop, 4*12))
        if n % 5 == 0:
            #plt.plot(generationList, AUPRCList, color = 'blue', label = 'AUPRC')
            plt.plot(generationList,averageMIList, color = 'orange', label = 'Average Percent MI of entire population')
            
            #Plots frequency of first tag
            plt.plot(generationList,idOrgFrequencyList[0], color = "purple", label = "Percent of population from" + id_list[0] + "organism")

            plt.xlabel("Generations of Simulation")
            plt.legend(loc = "upper left", bbox_to_anchor=(1,1))
            plt.ylim(0,1)
            plt.pause(0.000001)
        print(f"Generation {n}") 
    
    #Writes AUClist to file
    with open(output_file, mode = "w", newline = '') as output:
        csvwriter = csv.writer(output)
        csvwriter.writerow(["Generation", "Mean AUPRC", "Max AUPRC", "Average MI"] + [print(f"{tag} org frequency") for tag in id_list]) #Writes names of columns
        
        #Writes info to csv file 
        for i in range(len(generationList)):
            csvwriter.writerow([generationList[i],meanAUPRCList[i], maxAUPRCList[i], averageMIList[i]] + 
                               [idOrgFrequencyList[j][i] for j in len(idOrgFrequencyList)])
        output.close()
