# -*- coding: utf-8 

import json
import csv
import pickle
from scipy.stats import sem
from src.Objects.population import Population, fitMotifPopulation
from src.helperFunctions import reachedAsymptote

def updateLossDict(organismList, lossDict):
    """
    Helper function for simulation. Updates lossDict with statistics 
    about mean/stdev Information in the population

    """
    Average_IC = sum([obj.IC for obj in organismList])/len(organismList)
    SEM_IC = sem([obj.IC for obj in organismList])
    
    Average_MI = sum([obj.MI for obj in organismList])/len(organismList)
    SEM_MI = sem([obj.MI for obj in organismList])
    
    Average_TI = sum([obj.TI for obj in organismList])/len(organismList)
    SEM_TI = sem([obj.TI for obj in organismList])
    
    #Adds to lossDict
    lossDict["# Generation after mutation"].append(generationNum)
    lossDict["IC"].append(Average_IC); lossDict["SEM_IC"].append(SEM_IC)
    lossDict["MI"].append(Average_MI); lossDict["SEM_MI"].append(SEM_MI)
    lossDict["TI"].append(Average_TI); lossDict["SEM_TI"].append(SEM_TI)

if __name__ == "__main__":
    '''
    Json parameters:
        populationFile: str or null
          Path to pickle object of pregenerated population (null to generate the population)
          
        (Docstring of json parameters relavent for fitting the population is available in fitMotifPopulation function)
                            
        numGenerationsTerminate: int or null
            Number of generations to run information decay. Null will run simulation until an asymptotic value is reached
        
        outputCSV: str
            Path to output statistics about the population after each generation of mutation (information decay)
    '''
    
    #Opens the json file and reads
    with open("./src/settings_Figure1.json") as jfile:
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
    Eucledian_thres = json_file["Eucledian_Thres"]
    populationOutputFile = json_file["populationOutputFile"]
    populationDecayFile = json_file["populationDecayFile"]
    
    #Loads json input dealing with loss of information of population and num_generations in output file
    numGenerationsTerminate = json_file["NumGenerationTerminate"]
    outputCSV = json_file["lossOutputCSV"]
    
    
    if populationFile is None:
        fitMotifPopulation(num_organisms, num_BS, length_BS, gc_cont, mut_rate, read_in, MI_mode, target_TI_Column, Goal_IC, Goal_MI, Eucledian_thres, populationOutputFile)    
        
    else:
        #Else loads population file (after selection)
        with open(populationFile,'rb') as populationFileVar:
            population = pickle.load(populationFileVar)
    
    #Prepares list containing sublists with the elements in fields
    lossDict = {"# Generation after mutation":[], "IC":[], "MI":[], "TI":[], "SEM_IC":[], "SEM_MI":[], "SEM_TI":[]}

    fields = list(lossDict.keys())
    
    
    #Updates mutation rate
    for obj in population.population_list:
        obj.mut_rate = mut_rate
    
    #Counter variable for appending to lossList
    generationNum = 0

    
    #Runs until each generation reaches asymptote of TI
    if numGenerationsTerminate is None:
        while (not reachedAsymptote(lossDict["TI"])):
            #Calculates Average IC, MI, and TI of population and appends to lossDict
            updateLossDict(population.populationList, lossDict)
            
            #Mutates organism lines, updating Average IC, MI, and TI
            population.Lose_Information2()
            
            #Adds to generationNum
            generationNum += 1
            
    #Runs for a select number of generations
    else:
        while (generationNum <= numGenerationsTerminate):
            #Calculates Average IC, MI, and TI of population and appends to lossDict
            updateLossDict(population.populationList, lossDict)
            
            #Mutates organism lines, updating Average IC, MI, and TI
            population.Lose_Information2()
            
            #Adds to generationNum
            generationNum += 1
        
    
    #Prints to csv file
    with open(outputCSV, "w", newline = "") as outputcsv:
        csvwriter = csv.writer(outputcsv)
        
        #Write fields
        csvwriter.writerow(fields)
        
        #Writes list
        csvwriter.writerows(list(lossDict.values()))
          
    

    #Generates output variable if populationFile is not given
    if populationDecayFile is not None:
        with open(populationDecayFile, 'wb') as decayFile:
            pickle.dump(population, decayFile)
            decayFile.close()
        
        