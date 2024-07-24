'''
Script used to further evolve populations
'''
import pickle
import json

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
    Goal_IC = 12
    Goal_MI = 0.25
    
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
    
    with open(populationFile,'rb') as populationFileVar:
        population = pickle.load(populationFileVar)
    
    population.Goal_IC = Goal_IC
    population.Goal_MI = Goal_MI
    population.scrap_list = []
    population.MSE_list = []
    population.constrictMode = constrictMode
    
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
        
        print(population.MSE_list[-1])
        
    if populationOutputFile != "":
        with open(populationOutputFile, 'wb') as output:
            pickle.dump(population, output)
            output.close()