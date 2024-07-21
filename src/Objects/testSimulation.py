# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 22:40:21 2024

@author: thele
"""
import pickle

from src.Objects.coev_organism import coevPopulation
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
    testPop = coevPopulation(1000, 0.20, 0, 50, 2048, 32, 8, 32) 
    generationList = []
    meanAUPRCList = []
    MIList = []
    maxAUPRCList = []
    i = 0
    generationSum = {"Top AUPRC": 0}
    while generationSum["Top AUPRC"] < 1:
        #Dictionary with summary of statistics for generation
        generationSum = testPop.runOneGeneration()
        
        generationList.append(i)
        meanAUPRCList.append(generationSum["Mean AUPRC"])
        maxAUPRCList.append(generationSum["Top AUPRC"])
        
        MIList.append(generationSum["Mean Information"][1]/generationSum["Mean Information"][2])
        #avgHamList.append(averageHamming(testPop, 4*12))
        if i % 5 == 0:
            plt.plot(generationList, meanAUPRCList, color = 'blue', label = 'Average AUPRC')
            plt.plot(generationList,MIList, color = 'orange', label = 'Average Percent MI')
            plt.plot(generationList,maxAUPRCList, color = 'red', label = "Max AUPRC in population")
            plt.xlabel("Generations of Simulation")
            plt.legend()
            
            plt.pause(0.000001)
            
        
        print(
            f"Generation {i}: Mean Fitness - {generationSum['Mean Fitness']}, Mean % MI - {MIList[i]}, Mean AUPRC - {meanAUPRCList[i]}")


        i += 1
