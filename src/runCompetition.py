"""
Runs competition experiment with 1AUPRC organisms generated from fitCoevPopulation.py.

Loads in MI/IC coevOrganisms into a population, runs the coevolutionary simulation, 
and outputting the frequency of MI/IC organisms to a csv 

"""
import csv
import matplotlib.pyplot as plt
import json
import os
import sys
from src.Objects.coev_organism import coevPopulation, coevOrganism

def main():
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
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    #Grabs orgs
    input_orgs = []
    
    #Adds organisms to list with IDs
    for i in range(len(input_files)):
        importOrg = coevOrganism(importJSON = input_files[i])
        importOrg.orgTag = id_list[i]
        input_orgs.append(importOrg)
            
    
    motifPop = coevPopulation(popSize = num_org, targetMIproportion = 0.5, 
                              alpha = 0,beta=beta, importList = input_orgs)

    #Keeps track of results from each generation of coevolutionary simulation
    generationList = range(num_gen)
    #generationSummaryList = []
    freqOrgsList = {tag:[] for tag in id_list}
    
    
    #Runs simulation
    for n in range(num_gen):
        #Computes frequency statistics
        #Prepares dictionary that will contain counts of organisms with each ID
        freqOrgsID = {tag:0 for tag in id_list}
        
        #Iterates through populationList
        for i in range(motifPop.popSize):
            referenceOrg = motifPop.populationList[i]
            
            #Adds count to appropriate key
            for tag in freqOrgsID.keys():
                if referenceOrg.orgTag == tag:
                    freqOrgsID[tag] += 1
                    
        #Converts to frequency
        freqOrgsID = {tag:freqOrgsID[tag]/num_org for tag in id_list}
                   
        #Adds to list
        for tag in id_list:
            freqOrgsList[tag].append(freqOrgsID[tag])
            print(f"{tag}:{freqOrgsID[tag]}")
            
        #Goes through generation of evolutionary simulation
        motifPop.determineFitness()
        motifPop.Replicate()
        motifPop.Mutate()
        motifPop.Select()
        
        if n%25 == 0:
            print(f"Generation {n}")
        
        """if n%25 == 0:
            generationSum = motifPop.generationFitness()
            
            #Adds information to list
            generationList.append(n)
            generationSummaryList.append(generationSum)
            

            print(
                f"Generation {n}: Mean Fitness - {generationSum['Mean Fitness']}, Mean % MI - {generationSum['Mean BS %MI']}, Mean AUPRC - {generationSum['Mean AUPRC']}, Frequency of 70% MI - {freqOrgsID['MI_70']}")
            sys.stdout.flush()
        """
            
    
    #Writes AUClist to file
    with open(output_file, mode = "w", newline = '') as output:
        csvwriter = csv.writer(output)
        #rownames = list(generationSummaryList[0].keys())
        #csvwriter.writerow(rownames + [f"{tag} org frequency" for tag in id_list]) #Writes names of columns
        rownames = ["Generation"] + [f"{tag} org frequency" for tag in id_list]
        csvwriter.writerow(rownames)
        
        #Writes info to csv file 
        for i in range(len(generationList)):
            #csvwriter.writerow([generationSummaryList[i][rowname] 
                                #for rowname in rownames] + 
                               #[freqOrgsList[tag][i] for tag in id_list])
            csvwriter.writerow([generationList[i]] + [freqOrgsList[tag][i] for tag in id_list])
        output.close()
        
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
    main()
