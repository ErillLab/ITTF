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
    
    #Possible values: None, "IC", "MI"
    constrictMode = json_file["ConstrictMode"]
    if ((constrictMode != "IC") & (constrictMode != "MI")):
        constrictMode = None
        
        

    #Loads json input dealing with MSE threshold and output File creation
    Eucledian_thres = json_file["Eucledian_Thres"]
    populationOutputFile = json_file["populationOutputFile"]
    
    
    #Generates and fits motif population
    motifPop = MotifPopulation(num_organisms,num_BS,length_BS,gc_cont,mut_rate, Goal_IC, Goal_MI, constrictMode, MI_mode_org = MI_mode, read_in_org_list = read_in)
    
    aboveThreshold = True
    
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
                print(round(motifPop.EucledianDistance_list[i],4))

    

    #Exports to file path
    motifPop._export(populationOutputFile)

if __name__ == "__main__":
    '''
    Json parameters:
    ----------    
    num_organisms: int
        Number of motif organsims in the population
        
    num_BS: int
        Number of binding sites in the motif object
    
    length_BS: int
        Length of each binding site in the motif object
    
    gc_cont: float (0-1)
        GC content of motif. Weights random initialization and mutation 
    
    mut_rate: float (0-1)
        Expected proportion of all letters in the motif that are mutated in a generation
    
    read_in: list
        List of strings with paths to fasta files for initializing each MotifModel. 
        If not specified, defaults to none and randomly initialzes MotifModel
        
    MI_mode: int (1 or 2)
        Mode for computing MI of motifs. 1 = Overlapping pairs, 2 = Independent pairs
    
    target_TI_Column: float (0-2)
        Average number of bits of Total Information requested for each column
    
    Goal_IC: float (0-1)
        Proportion of information requested to be encoded through IC. 
        
    Goal_MI: float (0-1)
        Proportion of information requested to be encoded through MI. 1-Goal_IC 
         
    Eucledian_thres: float
        Eucledian distance threshold for fitness function.
        Must be reached by each member of the population to complete selection
        
    populationOutputFile: str or null
        Path to write fitted population to. 

    '''
    main()
    
    