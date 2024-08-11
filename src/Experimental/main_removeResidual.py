'''
Script used to further evolve populations
'''


import pickle
import json

def normalizeResiduals(populationList, MSE_thres):
    """
    Assumes first entry is MI
    Parameters
    ----------
    populationList : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if (populationList[0][-9:] != "MI.pickle"):
        raise ValueError("Put MI file first in list")
    
    with open(populationList[0], 'rb') as populationFileVar:
        populationMI = pickle.load(populationFileVar)
    
    #Computes average residual information
    residualInformationList = []
    for i in range(len(populationMI.population_list)):
        residualInformationList.append(populationMI.population_list[i].IC)
    goalResidual = sum(residualInformationList)/len(residualInformationList)
    
    populationIC = reduceResidual(populationList[1], goalResidual, MSE_thres)
    return populationIC
    
        
    

def reduceResidual(populationFile, Goal_MI, MSE_thres):
    #Opens the json file and reads 
    with open(populationFile,'rb') as populationFileVar:
        population = pickle.load(populationFileVar)
    
    population.Goal_MI = Goal_MI
    population.MSE_list = []
    population.scrap_list = []
    population.constrictMode = ""
    
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
        
    return population

    
if __name__ == "__main__":
    with open("./settings_residual.json") as jfile:
        json_file = json.load(jfile)
    
    populationList = json_file["populationList"]
    MSE_thres = json_file["MSE_thres"]
    
    editedPop = normalizeResiduals(populationList, MSE_thres)
    with open(populationList[1], 'wb') as output:
        pickle.dump(editedPop, output)
        output.close()

        
    


