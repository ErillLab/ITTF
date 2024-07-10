# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 22:40:21 2024

@author: thele
"""
import pickle

from coev_organism import coevPopulation
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
    testPop = coevPopulation(1000, 0.20, 0.50, 50, 2048, 4, 12, 32) 
    generationList = []
    AUPRCList = []
    MIList = []
    maxAUPRCList = []
    i = 0
    testPop.maxAUPRC = 0
    while testPop.maxAUPRC < 1:
        testPop.runOneGeneration()
        generationList.append(i)
        AUPRCList.append(testPop.meanAUPRC)
        maxAUPRCList.append(testPop.maxAUPRC)
        
        MIList.append(testPop.meanInformation[1]/testPop.meanInformation[2])
        #avgHamList.append(averageHamming(testPop, 4*12))
        if i % 5 == 0:
            plt.plot(generationList, AUPRCList, color = 'blue', label = 'AUPRC')
            plt.plot(generationList,MIList, color = 'orange', label = 'Average Percent MI')
            plt.plot(generationList,maxAUPRCList, color = 'red', label = "Max AUPRC in population")
            plt.xlabel("Generations of Simulation")
            plt.legend()
            
            plt.pause(0.000001)
            
        if i % 1000 == 0:
            with open(f"C:/Users/thele/anaconda3/envs/diPWM/src/Objects/populations/populationGeneration{i}.pickle", 'wb') as output:
                pickle.dump(testPop, output)
                output.close()
        
        print(
            f"Generation {i}: Mean Fitness - {testPop.meanFitness}, Mean % MI - {testPop.meanInformation[1]/testPop.meanInformation[2]}, Mean AUPRC - {testPop.meanAUPRC}")


        i += 1
