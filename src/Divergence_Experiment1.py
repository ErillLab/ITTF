# -*- coding: utf-8 
import json
import csv
import pickle
from src.Objects.population import Population
from scipy.stats import sem
from statistics import mean 
from copy import deepcopy
from src.helperFunctions import reachedAsymptote, AverageJSD
    

if __name__ == "__main__":
    '''
    Json parameters:
        populationFile: str 
          Path to pickle object of pregenerated population ("" to generate the population)
        
        num_organisms: int
            Number of motif organsims in the population

        target_TI: int ???
         Totalrepresenting the total TI the organisms should attempt to reach for
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
    
    
    if populationFile is None:
        
        #Initailizes population object with parameters
        population = Population(num_organisms,num_BS,length_BS,gc_cont,mut_rate, Goal_IC, Goal_MI, constrictMode = "", MI_mode_org = MI_mode, read_in_org_list = read_in)
        
        aboveThreshold = True
        
        #Repeats reproduce and select until metThreshold is true (when MSE values for each organism line is below threshold)
        while (aboveThreshold):
            
            #Sets threshold value to false (Will iterate through MSE values to find an above threshold case)
            aboveThreshold = False
            
            #Reproduces and selects 
            population.Reproduce()           
            population.Select()
            
            #Checks for a single organism with an MSE above the threshold
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
            
        
        