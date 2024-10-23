# -*- coding: utf-8 -*-
"""
Randomly initializes motifPopulation and fits organisms to encode MI/IC 
as specified in the settings_fitMotifPopulation JSON file

Used to generate populations for Information degradation, Divergence, 
and Genome Search experiments
"""
import json
from src.Objects.motifpopulation import MotifPopulation

def main():
    #Opens the json file and reads
    with open("./src/settings_fitMotifPop.json") as jfile:
        json_file = json.load(jfile)
        
    #Loads population (if specified)
    inputFile = json_file["inputFile"]
    
    if inputFile is None:
        #Loads json input into class variables
        num_organisms = json_file["Number_Organisms"]
        num_BS = json_file["Number_BindingSites"]
        length_BS = json_file["Length_BindingSites"]
        gc_cont = json_file["GC_cont"]
        mut_rate = json_file["Mutation_Rate"]
        read_in = json_file["ReadInFiles"]
    
    #Loads in parameters relavent to the fitness function
    MI_mode = json_file["MI_mode"]
    target_TI_Column = json_file["Target_TI_Column"]
    target_TI = target_TI_Column * length_BS
    Goal_IC = json_file["Prop_IC"] * target_TI
    Goal_MI = (1 - json_file["Prop_IC"]) * target_TI
    num_restartPop = json_file["restartPopulation"]
    
    #Possible values: None, "IC", "MI"
    constrictMode = json_file["ConstrictMode"]
    if ((constrictMode != "IC") & (constrictMode != "MI")):
        constrictMode = None
        
        

    #Loads json input dealing with MSE threshold and output File creation
    Eucledian_thres = json_file["Eucledian_Thres"]
    populationOutputFile = json_file["populationOutputFile"]
    
    
    #Generates motif population
    if inputFile is None:
        motifPop = MotifPopulation(num_organisms,num_BS,length_BS,gc_cont,
                                   mut_rate, Goal_IC, Goal_MI, constrictMode, 
                                   MI_mode_org = MI_mode, read_in_org_list = read_in)
    else:
        motifPop = MotifPopulation(importPath = inputFile)
        
        #Updates object variables
        motifPop.MI_mode_org = MI_mode
        motifPop.Goal_IC = Goal_IC
        motifPop.Goal_MI = Goal_MI
        motifPop.constrictMode = constrictMode
    
    
    #Fits motifPopulation
    aboveThreshold = True
    
    #Variables with information about population 
    completedOrganismList = []
    numOrg_fit = num_organisms
    generationNum = 0
    
    # Exit loop when ALL organisms have a MSE below threshold 
    #Repeats reproduce and select until metThreshold is true (when MSE values for each organism line is below threshold)
    while (aboveThreshold):
        
        #Sets threshold value to false (Will iterate through MSE values to find an above threshold case)
        aboveThreshold = False
        
        #Reproduces and selects 
        motifPop.Reproduce()           
        motifPop.Select()
        
        #Checks for a single organism with an ED above the threshold
        for i in range(len(motifPop.population_list)):
            if (motifPop.EucledianDistance_list[i] > Eucledian_thres):
                aboveThreshold = True
                numOrg_fit -= 1
                #print(round(motifPop.EucledianDistance_list[i],4))
            
                
        #Keeps track of generation num and number organisms that are below thresh
        generationNum += 1
        if generationNum % 100 == 0:
            print(f"Number Organisms below threshold: {numOrg_fit}")
            print(f"Generation num: {generationNum}")
                
        #Restarts population if option was selected and saves finished organisms
        if num_restartPop is not None:
            if generationNum % num_restartPop == 0:
                for i in range(len(motifPop.population_list)):
                    
                    #Grabs finished Organisms
                    if (motifPop.EucledianDistance_list[i] < Eucledian_thres):
                        completedOrganismList.append(motifPop.population_list[i].TFBS)
                
                #Checks if enough fitted organisms were found
                if len(completedOrganismList) >= num_organisms:
                    
                    #Reinitializes motifPop
                    motifPop = MotifPopulation(num_organisms,num_BS,length_BS,gc_cont,
                                               mut_rate, Goal_IC, Goal_MI, constrictMode, 
                                               MI_mode_org = MI_mode, 
                                               read_in_TFBS = completedOrganismList[:num_organisms])
                    
                    #Exits for loop
                    aboveThreshold = False
                    
                else:
                    print(f"Restarting population and saving {numOrg_fit} organisms")
                    #Restarts motifPop
                    motifPop = MotifPopulation(num_organisms,num_BS,length_BS,gc_cont,
                                               mut_rate, Goal_IC, Goal_MI, constrictMode, 
                                               MI_mode_org = MI_mode, read_in_org_list = read_in)    
        
        #Adjusts variable for next cycle
        numOrg_fit = num_organisms
            

    #Exports to file path
    motifPop._export(populationOutputFile)

if __name__ == "__main__":
    '''
    Json parameters:
    ----------    
    
    inputFile: null or str
        Set to null to randomly initialize population. Otherwise, loads in 
        motifPopulation from path (motifPopulation._export()) and runs simulation
        
    num_organisms: int or null (if inputfile specified)
        Number of motif organsims in the population
        
    num_BS: int or null (if inputfile specified)
        Number of binding sites in the motif object
    
    length_BS: int or null (if inputfile specified)
        Length of each binding site in the motif object
    
    gc_cont: float (0-1) or null (if inputfile specified)
        GC content of motif. Weights random initialization and mutation 
    
    mut_rate: float (0-1) or null (if inputfile specified)
        Expected proportion of all letters in the motif that are mutated in a generation
    
    read_in: list or null
        List of strings with paths to fasta files for initializing each MotifModel. 
        If not specified, defaults to none and randomly initialzes MotifModel
        
    MI_mode: int (1 or 2) 
        Mode for computing MI of motifs. 1 = Overlapping pairs, 2 = Independent pairs
    
    target_TI_Column: float (0-2)
        Average number of bits of Total Information requested for each column
    
    Prop_IC: float(0-1)
        Proportion of information to be encoded through IC. 
        Prop_MI = 1- PropMI
        
    ConstrictMode: str ("IC" or "MI") or null
        Fitness function of simulation either selects for solely "IC" or "MI"
        or both (null)
        
    Eucledian_thres: float
        Eucledian distance threshold for fitness function.
        Must be reached by each member of the population to complete selection
    
    restartPopulation: int or null
        If set to null, evolutionary simulation continues until 
        all organisms are below threshold. If set to int, runs 
        evolutionary simulation for # restartPopulation generations
        and saves the organisms that are below the threshold. 
        Restarts evolutionary simulation with new population and repeats
        until num_org organisms are generated below threshold.
        
        Helps get out of local minimas with specific organisms
        
    populationOutputFile: str 
        Path to write fitted population to. 

    '''
    main()
    
    