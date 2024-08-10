# -*- coding: utf-8 
import json
import csv
import pickle
from src.Objects.motifpopulation import fitMotifPopulation
from scipy.stats import sem
from statistics import mean 
from copy import deepcopy
from src.helperFunctions import reachedAsymptote, AverageJSD

def updateDictJSD(organismList, frozenMotifList, lossDict, generationNum):
    """
    Helper function for script. Computes average JSD of each organism
    and updates dict containing statistics about population
    
    ??? we talked about this in chat Elia, but I changed some of the functionality. 
    Let me know if I should still remove this function
    """
    
    #Depending on input in json file, adds to lossList
    JSDList = AverageJSD(frozenMotifList,[obj.TFBS for obj in organismList],2)
    
    #Updates lists
    lossDict["# Generation after mutation"].append(generationNum)
    lossDict["JSD"].append(mean(JSDList))
    lossDict["SEM_JSD"].append(sem(JSDList))
    
if __name__ == "__main__":
    '''
    Json parameters:
        populationFile: str or None
          Path to pickle object of pregenerated population (null to generate the population)
          
        (Docstring of json parameters relavent for fitting the population is available in fitMotifPopulation function)
                            
        numGenerationsTerminate: int or None
            Number of generations to run information decay. Null will run simulation until an asymptotic value is reached
        
        outputCSV: str
            Path to output statistics about the population after each generation of mutation (information decay)
    '''
    
    #Opens the json file and reads
    with open("./src/settings_divergence.json") as jfile:
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
    constrictMode = json_file["ConstrictMode"]
    
    #Loads in parameters relavent to the fitness function
    target_TI_Column = json_file["Target_TI_Column"]
    target_TI = target_TI_Column * length_BS
    Goal_IC = json_file["Prop_IC"] * target_TI
    Goal_MI = (1 - json_file["Prop_IC"]) * target_TI

    #Loads json input dealing with MSE threshold and output File creation
    Eucledian_thres = json_file["Eucledian_Thres"]
    populationOutputFile = json_file["populationOutputFile"]
    
    #Loads json input dealing with loss of information of population and num_generations in output file
    numGenerationsTerminate = json_file["NumGenerationTerminate"]
    outputCSV = json_file["lossOutputCSV"]
    
    
    if populationFile is None:
        fitMotifPopulation(num_organisms, num_BS, length_BS, gc_cont, mut_rate, read_in, MI_mode, target_TI_Column, Goal_IC, Goal_MI, constrictMode, Eucledian_thres, populationOutputFile)    
        
    else:
        #Else loads population file (after selection)
        with open(populationFile,'rb') as populationFileVar:
            population = pickle.load(populationFileVar)
    
    #Prepares dictionary containing sublists with the elements in fields
    lossDict = {"# Generation after mutation":[], "JSD":[],"SEM_JSD":[]}
    fields = list(lossDict.keys())
    frozenMotifList = deepcopy([obj.TFBS for obj in population.population_list])

    #Counter variable for appending to lossList
    generationNum = 0
    
    
    if numGenerationsTerminate is None:
        
        while (not reachedAsymptote(lossDict["JSD"])):   
            #Updates lists
            updateDictJSD([motifOrganism for motifOrganism in population.populationList], frozenMotifList, lossDict, generationNum)
            
            #Mutates organism lines, updating Average IC, MI, and TI
            population.Lose_Information2()
            
            #Adds to generationNum
            generationNum += 1
            
    #Runs for select number of gens
    else:
        while generationNum < numGenerationsTerminate:
            #Updates lists
            updateDictJSD([motifOrganism for motifOrganism in population.populationList], frozenMotifList, lossDict, generationNum)
            
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
        csvwriter.writerows(list(lossDict.values()))
        

    
            
    
    
        
        
    #Generates output variable if populationFile is not given
            
        
        