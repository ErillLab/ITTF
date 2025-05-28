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
from statistics import mean, stdev
import copy
import matplotlib.pyplot as plt
import json
import csv
import warnings
import sys
import itertools
import os

def splitList(inputList, numParts):
    a, b = divmod(len(inputList), numParts)
    return (inputList[i*a+min(i, b):(i+1)*a+min(i+1, b)] for i in range(numParts))

def finishedIteration(numGens,stopping_criteria, error_threshold, grid_threshold, topAUPRC, topFitness, meanAUPRC, meanMI, gridAUPRC, gridMI):
    """
    Helper function to determine when to end the simulation. 
    
    If stopping_criteria is an int (number generations required), checks to see if generation number has been reached
    If stopping_criteria is a str ("fitness") (Min L2 norm of fitness value from goal), checks to see if top fitness is close enough to goal
    If stopping_criteria is a str ("grid"), checks to see if population is on boundary of grid
    If stopping_criteria is null, checks to see if top AUPRC is 1
    
    Returns true if these conditions are met. False otherwise
    """
    if type(stopping_criteria) is int:
        return True if numGens >= stopping_criteria else False
    elif type(stopping_criteria) is str:
        if stopping_criteria == "fitness":
            return True if abs(topFitness) <= error_threshold else False
        elif stopping_criteria == "grid":
            #Computes boundaries
            topBoundary = gridAUPRC + grid_threshold
            bottomBoundary = gridAUPRC - grid_threshold
            leftBoundary = gridMI - grid_threshold
            rightBoundary = gridMI + grid_threshold
            
            if (meanAUPRC >= topBoundary) or (meanAUPRC <= bottomBoundary) or (meanMI <= leftBoundary) or (meanMI >= rightBoundary):
                return True 
            else: 
                return False
        else:
            raise Exception("If stopping_criteria is a string, must be fitness or grid")
    elif stopping_criteria is None:
        return True if topAUPRC == 1 else False
    else:
        raise Exception("Stopping Criteria must be an int, str, or None/null")
        
        
def computeSummaryStatistics(popDict,popSize):
    """
    Helper function to take split popDict returned from parallelized
    coevPopulation.determineFitness(), merge them, and compute summary statistics

    """
    #Computes summary statistics
    summaryDict = {}
    summaryDict["Mean Fitness"] = mean(popDict["Fitness List"])
    if inputPop.popSize > 1: 
        summaryDict["Stdev Fitness"] = stdev(popDict["Fitness List"])
    summaryDict["Top Fitness"] = max(popDict["Fitness List"])
    
    summaryDict["Mean AUPRC"] = mean(popDict["AUPRC List"])
    if inputPop.popSize > 1:
        summaryDict["Stdev AUPRC"] = stdev(popDict["AUPRC List"])
    summaryDict["Top AUPRC"] = max(popDict["AUPRC List"])
    
    
    #Variable for keeping track of mean Recognizer Information
    summaryDict["Mean Recognizer IC"] = mean(popDict["Recognizer IC"])
    summaryDict["Mean Recognizer MI"] = mean(popDict["Recognizer MI"])
    summaryDict["Mean Recognizer TI"] = mean(popDict["Recognizer TI"])
    summaryDict["Mean Recognizer %MI"] = mean(popDict["Recognizer %MI"]) 

    if popSize > 1:
        summaryDict["Stdev Recognizer IC"] = stdev(popDict["Recognizer IC"])
        summaryDict["Stdev Recognizer MI"] = stdev(popDict["Recognizer MI"])
        summaryDict["Stdev Recognizer TI"] = stdev(popDict["Recognizer TI"])
        summaryDict["Stdev Recognizer %MI"] = stdev(popDict["Recognizer %MI"])
    
    
    #Computes % Recognizer 
    #Variable for keeping track of mean binding site Information
    summaryDict["Mean BS IC"] = mean(popDict["BS IC"])
    summaryDict["Mean BS MI"] = mean(popDict["BS MI"])
    summaryDict["Mean BS TI"] = mean(popDict["BS TI"])
    summaryDict["Mean BS %MI"] = mean(popDict["BS %MI"])
    if popSize > 1:
        summaryDict["Stdev BS IC"] = stdev(popDict["BS IC"])
        summaryDict["Stdev BS MI"] = stdev(popDict["BS MI"])
        summaryDict["Stdev BS TI"] = stdev(popDict["BS TI"])
        summaryDict["Stdev BS %MI"] = stdev(popDict["BS %MI"])
    
    return summaryDict
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
        
        AUPRC: float
            Target AUPRC of population (top Organism)
        
        Beta: float
            Base of exponent for # of offspring. Affects selection vs drift of algorithm
        
        stop_criteria: int, str ("fitness" or "grid"), or None
            If int, number of generations to run simulation for. 
            If str, either stops after top organism reaches overcomes fitness threshold 
            or average AUPRC/MI reaches end of grid (If Grid_Threshold is specified)
            If none, runs until 1AUPRC organism is obtained. 
            
        error_threshold: float or none
            If stop_criteria is fitness, describes the max eucledian distance permitted
            by the top organism from the goal AUPRC/%MI before stopping the simulation
        
        Grid_Threshold: float or None
            Describes the distance from the end of the grid to the center (from AUPRC and targetMI)
            Only used if stop_criteria grid
            
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
    
    #AUPRC requested for recognizer
    req_AUPRC = json_file["AUPRC"]
    
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
    stopping_criteria = json_file["Stopping_Criteria"] 
    error_threshold = json_file["Error_Threshold"] 
    grid_threshold = json_file["Grid_Threshold"]
    run_mode = json_file["run_mode"]
    
    #Grabs output folder paths and names
    output_folder = json_file["Output_Folder"]
    output_csv_file_name = json_file["Output_CSV_FileName"]
    output_topAUPRC_org_file_name = json_file["Output_topAUPRC_Org_FileName"]
    output_topFitness_org_file_name = json_file["Output_topFitness_Org_FileName"]
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
            inputPop = coevPopulation(num_org, percent_MI, req_AUPRC, alpha, beta, genomeLength, numBS_Recognizer, lenBS, numBS,importList = [coevOrganism(importJSON = filePath) for filePath in input_files])
        
        #Initializes population from file path
        elif type(input_files) is str:
            inputPop = coevPopulation(importFile = input_files)
        
        #Random Initialization
        else:
            inputPop = coevPopulation(num_org, percent_MI, req_AUPRC, alpha, beta, genomeLength, numBS_Recognizer, lenBS, numBS) 
       
        #Computes fitness for first iteration
        popDict = inputPop.determineFitness()
        generationSum = computeSummaryStatistics(popDict,inputPop.popSize)
        generationSummaryList.append(generationSum)
        
        #Grabs grid's AUPRC/MI value for grid experiment (not used if not running grid experiment)
        gridAUPRC = generationSum["Mean AUPRC"]
        gridMI = generationSum["Mean BS %MI"]
        
    else:
        inputPop = None       

    finishedIterating = False
    
    #Runs parallelized iterations of coevolutionary simulation
    while not finishedIterating:        
        
        if rank ==0:
            #Replicates
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
        
        #Computes fitness 
        if rank == 0:
            #Splits population among nodes
            data = splitList(inputPop.populationList,size)
            
        else: 
            data = None
        
        #Initializes subpopulations and computes fitness
        local_data = comm.scatter(data, 0)
        input_subpop = coevPopulation(len(local_data), percent_MI, req_AUPRC,alpha, beta, 
                                      genomeLength, numBS_Recognizer, lenBS, 
                                      numBS,importList = local_data)
        subpopDict = input_subpop.determineFitness() #Grabs lists about organisms
        
        #Gathers whole population in root note
        gathered_data = comm.gather((input_subpop,subpopDict),root=0)
        if rank == 0:
            orgList = []
            fitnessList = []
            popDictList = []
            for item in gathered_data:
                motifPop = item[0]
                orgList += motifPop.populationList
                fitnessList += motifPop.fitnessList
                popDictList.append(item[1])
            
            #Change name of object to uppercase (to indicate it's a class instead of function)
            inputPop = coevPopulation(len(inputPop.populationList), percent_MI, req_AUPRC,alpha, beta, 
                                          genomeLength, numBS_Recognizer, lenBS, 
                                          numBS,importList = orgList)
            inputPop.fitnessList = fitnessList
            
            #Merges subdictionaries
            popDict = popDictList[0]
            for i in range(1,len(popDictList)):
                for key in list(popDict.keys()):
                    popDict[key] += popDictList[i][key]
            
            
            #Computes summary statistics
            generationSum = computeSummaryStatistics(popDict, inputPop.popSize)
            generationSummaryList.append(generationSum)
            
            #Grabs top fitness, topAUPRC, meanAUPRC, and meanMI (stopping criteria)
            topFitness = generationSum["Top Fitness"]
            topAUPRC = generationSum["Top AUPRC"]
            meanAUPRC = generationSum["Mean AUPRC"]
            meanMI = generationSum["Mean BS %MI"]
            
            #Grabs top performing organisms (stopping criterion)
            topAUPRCindex = popDict["AUPRC List"].index(topAUPRC)
            topFitnessindex = popDict["Fitness List"].index(topFitness)
            
            
            
            numGens += 1

            if numGens%25 == 0:
                print(
                    f"Generation {numGens}: Mean Fitness - {generationSum['Mean Fitness']}, Mean % MI - {generationSum['Mean BS %MI']}, Mean AUPRC - {generationSum['Mean AUPRC']}")
                sys.stdout.flush()
            
            
            finishedIterating = finishedIteration(numGens,stopping_criteria, error_threshold, grid_threshold, topAUPRC, topFitness, meanAUPRC, meanMI,gridAUPRC,gridMI)
            print(numGens,stopping_criteria, error_threshold, topAUPRC, topFitness)
        
        finishedIterating = comm.scatter([finishedIterating for i in range(size)], root=0)                      
        
    if rank ==0:
        #Defines full paths unless provided output is absolute path
        output_csv_file = output_csv_file_name if os.path.isabs(output_csv_file_name) else output_subfolder + output_csv_file_name
        output_topAUPRC_org_file = output_topAUPRC_org_file_name if os.path.isabs(output_topAUPRC_org_file_name) else output_subfolder + output_topAUPRC_org_file_name
        output_topFitness_org_file = output_topFitness_org_file_name if os.path.isabs(output_topFitness_org_file_name) else output_subfolder + output_topFitness_org_file_name
        output_pop_file = output_pop_file_name if os.path.isabs(output_pop_file_name) else output_subfolder + output_pop_file_name
        csv_exists = os.path.exists(output_csv_file)
        
        # Check if CSV has only headers or is empty
        write_to_csv = True
        if csv_exists:
            with open(output_csv_file, mode="r", newline="") as f:
                reader = list(csv.reader(f))
                if len(reader) > 1:  # More than just header
                    write_to_csv = False
                    warnings.warn("CSV file already written to.")
            
        #Saves CSV file
        if output_csv_file and write_to_csv:
            with open(output_csv_file, mode = "a", newline = "") as output:
                csvwriter = csv.writer(output)
                
                rownames = list(generationSummaryList[0].keys())
                if not csv_exists:
                    csvwriter.writerow(rownames)
               
                #Writes info to csv file 
                for i in range(len(generationSummaryList)):
                    csvwriter.writerow([generationSummaryList[i][rowname] 
                                        for rowname in rownames])
                output.close()

            
        #Grabs Best classifier and most fit organism and saves to file
        topAUPRCorg = inputPop.populationList[topAUPRCindex]
        topFitnessorg = inputPop.populationList[topFitnessindex]
        
        if output_topAUPRC_org_file:
            topAUPRCorg._export(output_topAUPRC_org_file)
            
        if output_topFitness_org_file:
            topFitnessorg._export(output_topFitness_org_file)            
        
        #Saves population file
        if output_pop_file:
            inputPop._export(output_pop_file)
        
    
    else:
        print(f"[Process {rank} completed]")

   