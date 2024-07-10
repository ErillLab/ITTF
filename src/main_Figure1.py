# -*- coding: utf-8 

import json
import csv
import pickle
from scipy.stats import sem
from Objects.population import Population
from helperFunctions import reachedAsymptote


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
    with open("./settings_Figure1.json") as jfile:
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
    populationDecayFile = json_file["populationDecayFile"]
    
    #Loads json input dealing with loss of information of population and num_generations in output file
    numGenerationsPrint = json_file["NumGenerationPrint"]
    numGenerationsTerminate = json_file["NumGenerationTerminate"]
    outputCSV = json_file["lossOutputCSV"]
    
    
    if populationFile == "":
        
        #Initailizes population object with parameters
        population = Population(num_organisms,num_BS,length_BS,gc_cont,mut_rate, Goal_IC, Goal_MI, constrictMode, MI_mode_org = MI_mode, read_in_org_list = read_in)
        
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
            
            print("Still running:", population.MSE_list[-1])
        
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
    fields = ["# Generation after mutation", "IC", "MI", "TI", "SEM_IC", "SEM_MI", "SEM_TI"]
    
    #Calculates Average IC, MI, and TI of population
    Average_IC = sum([obj.IC for obj in population.population_list])/len(population.population_list)
    SEM_IC = sem([obj.IC for obj in population.population_list])
    
    Average_MI = sum([obj.MI for obj in population.population_list])/len(population.population_list)
    SEM_MI = sem([obj.MI for obj in population.population_list])
    
    Average_TI = sum([obj.TI for obj in population.population_list])/len(population.population_list)
    SEM_TI = sem([obj.TI for obj in population.population_list])
    
    #Updates mutation rate
    for obj in population.population_list:
        obj.mut_rate = mut_rate
    
    #Counter variable for appending to lossList
    generationNum = 0

    #Runs until each generation reaches asymptote of TI
    #while (generationNum <= numGenerationsTerminate):
    while (not reachedAsymptote(abrevLossList)):
        
        #Depending on input in json file, adds to lossList
        if ((generationNum % numGenerationsPrint) == 0):
            lossList.append([generationNum, Average_IC,Average_MI, Average_TI, SEM_IC, SEM_MI, SEM_TI])
            abrevLossList.append(Average_TI)
            
        #Mutates organism lines, updating Average IC, MI, and TI
        population.Lose_Information2()
        #print(population.population_list[0].IC)
        Average_IC = sum([obj.IC for obj in population.population_list])/len(population.population_list)
        SEM_IC = sem([obj.IC for obj in population.population_list])
        
        Average_MI = sum([obj.MI for obj in population.population_list])/len(population.population_list)
        SEM_MI = sem([obj.MI for obj in population.population_list])
        
        Average_TI = sum([obj.TI for obj in population.population_list])/len(population.population_list)
        SEM_TI = sem([obj.TI for obj in population.population_list])    
        
        #Adds to generationNum
        generationNum += 1
        
    
    #Prints to csv file
    with open(outputCSV, "w", newline = "") as outputcsv:
        csvwriter = csv.writer(outputcsv)
        
        #Write fields
        csvwriter.writerow(fields)
        
        #Writes list
        csvwriter.writerows(lossList)
          
    

    #Generates output variable if populationFile is not given
    if populationDecayFile != "":
        with open(populationDecayFile, 'wb') as decayFile:
            pickle.dump(population, decayFile)
            decayFile.close()
        
        