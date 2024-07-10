# -*- coding: utf-8 -*-
import pickle
import csv
import os
from scipy.stats import sem
from statistics import mean

def reachedAsymptote(yList):
    """
    Accepts 2D lists where each successive element (list) is the next generation of mutations
    
    yxList - 1D List where the first element of each list is the x value and the second element is the y value
                                                                  
    Called by loop in the main function used to induce mutation in the population. Returns true when slope hits threshold close to 0 (Asymptote).
    Returns false is slope does not hit threshold or 10 generations have not yet been generated
                                                                  
    """
    
    #Edge case
    if (len(yList) <= 30):
        return False
    
    #Computes Average of last 10 elements of list and the 10 before it
    averageYValcurr = mean([x for x in yList[-10:]])    
    averageYvalprev = mean([x for x in yList[-30:-19]])
    
    #Computes the magnitude of the difference
    differenceMagnitude = abs(averageYValcurr - averageYvalprev)
    print(differenceMagnitude)
    #Returns true if reached asymptote (slope close to 0)
    if (differenceMagnitude <= 0.0005):
        #Continues iteration until number of generations is a multiple of 100
        if ( ((len(yList) - 1)%100) == 0):
            return True
        else:
            return False
    else:
        return False
    

def equalizeGenerations(decayFile, csvFile, numGenerations):
    """
    Helper function to take a decayFile associated with the numbers in the csvFile and further decay the population till numGenerations

    Parameters
    ----------
    decayFile : string
        pickle file path containing decayed population
    csvFile : string
        csv file path containing metrics from pickle population in decayFile
    numGenerations : int
        number of generations to decay till

    Returns
    -------
    None.

    """
    
    #Loads population
    with open(decayFile,'rb') as populationFileVar:
        population = pickle.load(populationFileVar)
        
        populationFileVar.close()
    
    #Recreates lossList
    lossList = []
    
    #Loads csv file
    with open(csvFile, "r", newline = "") as csvFileHandle:
        csvreader = csv.reader(csvFileHandle)
    
        for row in csvreader:
            lossList.append(row)
            
        csvFileHandle.close()
        
    #Gets last generation number
    latestGeneration = int(lossList[-1][0])
    
    #Induces mutation until we get to same generation number as numGenerations
    while latestGeneration < numGenerations:
        
        #Mutates organism lines, updating Average IC, MI, and TI
        population.Lose_Information()
        Average_IC = sum([obj.IC for obj in population.population_list])/len(population.population_list)
        SEM_IC = sem([obj.IC for obj in population.population_list])
        
        Average_MI = sum([obj.MI for obj in population.population_list])/len(population.population_list)
        SEM_MI = sem([obj.MI for obj in population.population_list])
        
        Average_TI = sum([obj.TI for obj in population.population_list])/len(population.population_list)
        
        #Adds to generation counter
        latestGeneration += 1
        
        #Adds to list
        lossList.append([latestGeneration, Average_IC, Average_MI, Average_TI, SEM_IC, SEM_MI])
        
    
    #After decay, resaves csv file
    with open(csvFile, "w", newline = "") as outputcsv:
        csvwriter = csv.writer(outputcsv)
                
        #Writes list
        csvwriter.writerows(lossList)
        
        outputcsv.close()
    
def equalizeFiles(filePairs):
    """
    Function passed a list of file pairs (1st = decay pickle population file and 2nd = CSVfile). 
    This function determines the maxmium generation among all of the csv files and induces additional 
    mutations so that all csv files have the same # of generations

    Parameters
    ----------
    filepairs: List of tuples
        Collection of tuples where each tuple contains #1 the decayed population file path and #2 the csv file path
        

    Returns
    -------
    None.

    """    
    
    #Stores csv files in the following list
    csvFileList = []
    
    #Variable for max generation number
    maxGeneration = 0
    
    #Populates list    
    for i in range(len(filePairs)):
        csvFileList.append(filePairs[i][1])
    
    #Iterates through list looking for largest num of generations in csv file
    for filePath in csvFileList:
        with open(filePath, "r", newline = "") as csvFileHandle:
            csvreader = csv.reader(csvFileHandle)
            
            csvGeneration = int(list(csvreader)[-1][0])
            #Grabs last row and checks if the value is greater than previous maxGeneration
            if (csvGeneration) > maxGeneration:
                maxGeneration = csvGeneration
            csvFileHandle.close()
            
    #Runs equalizeGenerations on each tuple of files
    for i in range(len(filePairs)):
        equalizeGenerations(filePairs[i][0], filePairs[i][1], maxGeneration)
        
def equalizeSubFiles(DecayPopulationDir, csvFolderDir):
    """
    Runs equalizeFiles on all populations contained 

    Parameters
    ----------
    csvFolderDir : TYPE
        DESCRIPTION.
    DecayPopulationDir : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    filesList = []
    for subdir, dirs, files in os.walk(csvFolderDir):
        for file in files:
            noExtensionFile = os.path.splitext(file)[0]
            filesList.append((DecayPopulationDir + noExtensionFile + ".pickle", csvFolderDir + noExtensionFile + ".csv"))
    
    equalizeFiles(filesList)
            
