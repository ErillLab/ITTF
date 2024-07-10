import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import pickle
from coev_organismReal import coevPopulation
import matplotlib.pyplot as plt
from scipy.spatial.distance import hamming
from statistics import mean

def averageHamming(coevPop, geneLength):
    geneList = []
    HammingList = []
    for i in coevPop.populationList:
        geneList.append(i.genomeObj.genome[:geneLength])
    for i in range(len(geneList)-1):
        tempList = []
        
        for j in range(i+1,len(geneList)):
            tempList.append(hamming(list(geneList[i]), list(geneList[j])))
        HammingList.append(mean(tempList))
    
    return mean(HammingList) 
if __name__ == "__main__":
    with open("C:/Users/thele/Research Projects/testRepo/src/Objects/populations/populationpseudoExp3_v3.pickle", 'rb') as inputFile:
        inputPop = pickle.load(inputFile)
        inputFile.close()

    generationList = []
    AUPRCList = []
    MIList = []
    maxAUPRCList = []
    MItagProportionList = []
    MI80percentMI = []
    MI20percentMI = []
    
    inputTags = ["MI_80", "MI_20"]
    
    inputPop.maxAUPRC = 0
    for n in range(10000):
        #Runs one generation
        inputPop.runOneGeneration()
        
        #Adds to lists for mean AUPRC, max AUPRC, and mean information
        generationList.append(n)
        AUPRCList.append(inputPop.meanAUPRC)
        maxAUPRCList.append(inputPop.maxAUPRC)
        MIList.append(inputPop.meanInformation[1]/inputPop.meanInformation[2])
                   
        #Prepares list containing information about each ID in population
        idLists = [[] for i in range(len(inputTags))]
            
        for i in range(inputPop.popSize):
            referenceOrg = inputPop.populationList[i]
            for k in range(len(inputTags)):
                if referenceOrg.orgID == inputTags[k]:
                    idLists[k].append(referenceOrg)
        
        averageMIforgroup = [] #Note this is the average % MI in each tag subpop
        #2D list
        for elementList in idLists: 
            if (len(elementList) == 0):
                averageMIforgroup.append(0)
            else:
                averageMIforgroup.append(mean([(org.organismInformation()[1])/(org.organismInformation()[2]) for org in elementList]))
            
        #Adds MI_80 tag to list
        MItagProportionList.append(len(idLists[0])/(1000))
        #Adds % MI for each subpopulation
        MI80percentMI.append(averageMIforgroup[0])
        MI20percentMI.append(averageMIforgroup[1])
        
                
        
        
        #avgHamList.append(averageHamming(testPop, 4*12))
        if n % 5 == 0:
            #plt.plot(generationList, AUPRCList, color = 'blue', label = 'AUPRC')
            plt.plot(generationList,MIList, color = 'orange', label = 'Average Percent MI of entire population')
            #plt.plot(generationList,maxAUPRCList, color = 'red', label = "Max AUPRC in population")
            plt.plot(generationList, MItagProportionList, color = "purple", label = "Percent of population from MI_80 organism")
            plt.plot(generationList, MI80percentMI, color = "brown", label = "Percent MI for MI_80 subpopulation")
            plt.plot(generationList, MI20percentMI, color = "black", label = "Percent MI for MI_20 subpopulation")

            plt.xlabel("Generations of Simulation")
            plt.legend(loc = "upper left", bbox_to_anchor=(1,1))
            plt.ylim(0,1)
            plt.pause(0.000001)
        print(f"Generation {n}")  
