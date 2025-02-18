# -*- coding: utf-8 -*-
"""
Script used to fit an inputted coevPopulation (random initialization, passed in by file, etc.).
Runs the coevolutionary simulation, selecting for offspring with recognizers 
that better classify their binding sites

Used to:
    - Generate a 1AUPRC coevOrganism to run the competition experiment (runCompetition.py)
    - Drift coevPopulations for the drift experiment 
"""
from src.Objects.coev_organism import coevOrganism, coevPopulation
from src.helperFunctions import check_dir, def_out_dir_path
import copy
import matplotlib.pyplot as plt
import json
import csv
import time
import sys
import itertools

def splitList(inputList, numParts):
    a, b = divmod(len(inputList), numParts)
    return (inputList[i*a+min(i, b):(i+1)*a+min(i+1, b)] for i in range(numParts))

def finishedIteration(genNum, numTerminate, topAUPRC):
    """
    Helper function to determine when to end the simulation. If numTerminate is specified,
    returns True if genNum > numTerminate. Otherwise, returns True if topAUPRC ==1
    """
    if numTerminate:
        return True if genNum > numTerminate else False
    else:
        return True if topAUPRC >= 1 else False
    

if __name__ == "__main__":
    """
    JSON_File parameters:
        Input_Files: list, str, or None
            Set to list to initialize population from list of coev_organism paths. Set to str 
            to initialize population from file path. Set to none to randomly initialize 
            population
        
        Number_Organisms: int
            Number of organisms in coevolutionary simulation population. Duplicates
            organisms in Input_files.
        
        ... (Parameters for coevPopulation. For a more detailed overview, take a look at src/Objects/coev_organisms)
        
        Beta: float
            Base of exponent for # of offspring. Affects selection vs drift of algorithm
        
        Number_Generations: int or None
            Number of generations to run simulation for. 
            If none, runs until 1AUPRC organism is obtained. 
        
        
        run_mode: str
            "serial" or "parallel". If Parallelized, saves results to subfolders
            in directory written in Output_Folder
            
        Output_Folder: str 
            Path to output runs to. Runs will be stored in subfolders in
            this directory
            
        Output_CSV_FileName: str or None
            File name for statisticsc about the population after
            each generation of simulation (Max AUPRC, Mean AUPRC, Mean MI). 
            Saved in subdirectory in Output_Folder. Set to none to not save 
            
        Ouput_Org_FileName: str or None
            File name to save 1AUPRC organism to (JSON). If 1AUPRC organism has not been found, 
            outputs top organism. Set Num_Generations to none to guarantee 1AUPRC organism is found.
            Saved in subdirectory in Output_Folder. Set to none to not save 

        
        Output_Pop_FileName: str or None
            File Name for output population (JSON). 
            Saved in subdirectory in Output_Folder. Set to none to not save
            
    """
    #Temporary to limit time of run
    start_time = time.time()
    max_time_hours = 97
    max_time = max_time_hours * 60 * 60
    
    #Opens the json file and reads
    with open("./src/settings_fitCoevPop.json") as jfile:
        json_file = json.load(jfile)
    
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    #Input File to run simulation on pregenerated population (Used to resume coevolutionary simulation or change parameters. None for random initialization)
    input_files = json_file["Input_Files"]
    
    #Number of organisms for coevolutionary population
    num_org = json_file["Number_Organisms"]
    
    #Percent MI requested for recognizer
    percent_MI = json_file["Percent_MI"]
    
    #Parameters for fitness functions
    #Alpha - Weighting on Percent MI (1 = Equal weight as AUPRC)
    #Beta - Base of exponent for computing num offspring
    alpha = json_file["Alpha"]
    beta = json_file["Beta"]
    
    #Settings for each coevolutionary organism. Must only be specified if input_files is None
    genomeLength = json_file["Genome_Length"]
    numBS_Recognizer = json_file["Number_BindingSites_Recognizer"]
    lenBS = json_file["Length_BindingSites"]
    numBS = json_file["Number_BindingSites"]
    
    #Grabs run parameters
    num_gen = json_file["Number_Generations"]
    run_mode = json_file["run_mode"]
    
    #Grabs output folder paths and names
    output_folder = json_file["Output_Folder"]
    output_csv_file_name = json_file["Output_CSV_FileName"]
    output_org_file_name = json_file["Output_Org_FileName"]
    output_pop_file_name = json_file["Output_Pop_FileName"]
    output_diversity_file_name = json_file["Output_Diversity_FileName"]
    
    #Definites output path with run mode
    if rank == 0:
        output_subfolder = def_out_dir_path(run_mode, output_folder)
        check_dir(output_subfolder)
    
    #Statistics of population throughout simulation
    generationList = []
    generationSummaryList = []
    meanAUPRCList = []
    sdAUPRCList = []
    maxAUPRCList = []
    
    #Recognizer information
    averageMIList_recognizer = []
    sdMIList_recognizer = []
    averageTIList_recognizer = []
    sdTIList_recognizer = []
    
    #Binding site information
    averageMIList_bs = []
    sdMIList_bs = []
    averageTIList_bs = []
    sdTIList_bs = []    
    
    
    #Counter variables
    numGens=0
    print(rank,":",size)
    
    if rank == 0:
        #Initializes population from list of organisms
        if type(input_files) is list:
            inputPop = coevPopulation(num_org, percent_MI, alpha, beta, genomeLength, numBS_Recognizer, lenBS, numBS,importList = [coevOrganism(importJSON = filePath) for filePath in input_files])
        
        #Initializes population from file path
        elif type(input_files) is str:
            inputPop = coevPopulation(importFile = input_files)
        
        #Random Initialization
        else:
            inputPop = coevPopulation(num_org, percent_MI, alpha, beta, genomeLength, numBS_Recognizer, lenBS, numBS) 
       
    
    else:
        inputPop = None
        
    #Initialized for first iteration
    generationSum = {"Top AUPRC": 0}
            
    #Runs parallelized iterations of coevolutionary simulation
    while not finishedIteration(numGens,num_gen,generationSum["Top AUPRC"]):        
        if rank == 0:
            #Splits population among nodes
            data = splitList(inputPop.populationList,size)
            
        else: 
            data = None
        
        #Initializes subpopulations and computes fitness
        local_data = comm.scatter(data, 0)
        input_subpop = coevPopulation(len(local_data), percent_MI, alpha, beta, 
                                      genomeLength, numBS_Recognizer, lenBS, 
                                      numBS,importList = local_data)
        input_subpop.determineFitness()
        
        #Gathers whole population in root note
        gathered_data = comm.gather(input_subpop)
        if rank == 0:
            orgList = []
            fitnessList = []
            for motifPop in gathered_data:
                orgList += motifPop.populationList
                fitnessList += motifPop.fitnessList
            
            inputPop = coevPopulation(len(inputPop.populationList), percent_MI, alpha, beta, 
                                          genomeLength, numBS_Recognizer, lenBS, 
                                          numBS,importList = orgList)
            inputPop.fitnessList = fitnessList
            #print(len(fitnessList), f"Length of gathered_data: {len(gathered_data)}, and 2nd length: {len(gathered_data[0])}")
            inputPop.Replicate() 
        
        #Sends out future offspring to nodes to induce mutations
        if rank == 0:
            data_2 = splitList(inputPop.replicatedList,size)
        else:
            data_2 = None
            
        local_data_2 = comm.scatter(data_2,0)
        for k in range(len(local_data_2)):
            local_data_2[k].induceMutation()
        gathered_data_2 = comm.gather(local_data_2)       
        
        #Selects offspring for next generation
        if rank==0:
            inputPop.replicatedList = list(itertools.chain(*gathered_data_2))
            #Mutates population
            inputPop.Select()
            
            #Saves population statistics 
            if numGens%25 == 0:
                numCycles = numGens // 25
                generationSum = inputPop.generationFitness()
                diversityList = copy.deepcopy(generationSum["Parent organisms"])
                del generationSum["Parent organisms"]
                
                #Adds information to list
                generationList.append(numGens)
                generationSummaryList.append(generationSum)
                
                #Saves to csv file
                with open(output_subfolder + output_diversity_file_name, 'a', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(diversityList)
                """
                meanAUPRCList.append(generationSum["Mean AUPRC"])
                sdAUPRCList.append(generationSum["Stdev AUPRC"])
                maxAUPRCList.append(generationSum["Top AUPRC"])
                
                #Grabs recognizer Information
                mean_IC, mean_MI, mean_TI = generationSum["Mean Recognizer Information"]
                averageMIList.append(mean_MI/mean_TI)
                """
                print(
                    f"Generation {numGens}: Mean Fitness - {generationSum['Mean Fitness']}, Mean % MI - {generationSum['Mean BS %MI']}, Mean AUPRC - {generationSum['Mean AUPRC']}")
                sys.stdout.flush()
                
        numGens += 1
            
    
    if rank ==0:
        #Defines full paths
        output_csv_file = output_subfolder + output_csv_file_name
        output_org_file = output_subfolder + output_org_file_name
        output_pop_file = output_subfolder + output_pop_file_name
        
        #Saves CSV file
        if output_csv_file:
            with open(output_csv_file, mode = "w", newline = '') as output:
                csvwriter = csv.writer(output)
                rownames = list(generationSummaryList[0].keys())
                csvwriter.writerow(rownames) #Writes names of columns
                
                #Writes info to csv file 
                for i in range(len(generationList)):
                    csvwriter.writerow([generationSummaryList[i][rowname] 
                                        for rowname in rownames])
                output.close()
            
        #Grabs Best classifier and saves to file 
        for i in range(len(inputPop.populationList)):
            if (inputPop.populationList[i].returnAUPRC() == generationSum["Top AUPRC"]):
                topOrganism = inputPop.populationList[i]
        
        if output_org_file:
            topOrganism._export(output_org_file)
        
        #Saves population file
        if output_pop_file:
            inputPop._export(output_pop_file)

   