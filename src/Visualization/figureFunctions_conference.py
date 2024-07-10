# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 22:14:46 2024

@author: thele
"""
import os
os.chdir("/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src")
from sklearn.metrics import auc
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import log2, log, e

Figure1Files = ["./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_1IC.csv","./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_1MI.csv" ]

SuppFigure1Files = ["./Output/Loss/Paper/SuppFigure_1/population_100_128_2_1bit_0IC.csv", "./Output/Loss/Paper/SuppFigure_1/population_100_128_2_1bit_0.25IC.csv",
            "./Output/Loss/Paper/SuppFigure_1/population_100_128_2_1bit_0.5IC.csv", "./Output/Loss/Paper/SuppFigure_1/population_100_128_2_1bit_0.75IC.csv",
            "./Output/Loss/Paper/SuppFigure_1/population_100_128_2_1bit_1IC.csv"]

Figure2Files = [
    ["./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_0.5IC.csv","./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_0.5MI.csv"],
    ["./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_1IC.csv", "./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_1MI.csv"],
    ["./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_1.5IC.csv", "./Output/Loss/Paper/Figure_2/population_100_128_2_1bit_1.5MI.csv"]]

Figure3Files = [
    "./Output/Loss/Paper/Figure_3/population_100_128_12_0.5_1IC.csv",
    "./Output/Loss/Paper/Figure_3/population_100_128_12_0.5_1MI.csv"
    ]

SuppFigure2Files = [
    "./Output/Loss/Paper/SuppFigure_2/population_100_128_2_1bit_1IC.csv",
    "./Output/Loss/Paper/SuppFigure_2/population_100_128_2_1bit_1MI.csv"
    ]

Figure4Files = [
    "./Output/Loss/GenomeExperiment_16bp_12bit_SmallGenome/population_5_128_16_0.75_1IC.csv",
    "./Output/Loss/GenomeExperiment_16bp_12bit_SmallGenome/population_5_128_16_0.75_1MI.csv"
    ]


        
def computeAUC(xVec, yVec, slidingWindow):
    """
    Helper function for figure generation. 
    Computes AUC via the trapezoidal rule. Removes asymptotic area under the curve. 

    xVec: list - list of elements on xaxis. Should match the length of yVec
    yVec: list - list of elements on yaxis. Shuld match the length of xVec
    slidingWindow: int - Length of sliding window for computing AUC

    Returns
    -------
    Area under the curve
    """
    
    #Edge case where dataset isn't big enough for sliding window
    if (len(yVec) <= slidingWindow):
        raise ValueError("Not enough generations of mutations to compute AUC")
    
    #Initializes variable holding asymptotic TI value
    asymptoteTI = -1
    
    #Finds asymptote 
    for i in range(len(yVec) - slidingWindow):
        
        #Grabs values for slope calculation
        x1 = xVec[i]; x2 = xVec[i+(slidingWindow - 1)]
        y1 = yVec[i]; y2 = yVec[i+(slidingWindow-1)]
        
        slope = abs((y2-y1)/(x2-x1))
        
        if slope <= 0.00005:
            asymptoteTI = yVec[i]
    
    #If asymptote never found, assumes last TI value is asymptote
    if asymptoteTI == -1:
        asymptoteTI = yVec[-1]
    
    #Computes AUC using scikit
    AUC = auc(np.array(xVec), np.array(yVec))
    
    #computes area under asymptote 
    AUCasymptote = asymptoteTI * len(xVec)
    
    #Returns AUC below information curve but above asymptote
    return (AUC - AUCasymptote)

def plotGenericInformationCurve(fileList, text = True, pltTitle = None, xLabel = None, ylabel = None, legendLabel = None):
    """
    Plotting function returning the seaborn figureplot
    Plots MI and IC curves of two different populations. Called by several figure functions
    
    filelist - 2 element list containing strings describing relative path to IC and MI csv files
    text - Boolean indicating whether to include text in the graph
    Returns
    -------
    figureplot

    """
    #Populates dataframe from csv files and creates combined data frame
    ICdf = pd.read_csv(fileList[0])
    MIdf = pd.read_csv(fileList[1])

    #Manipulates combined dataframe for plotting purposes
    combineddf = pd.DataFrame({
        'Generation': np.concatenate((ICdf["# Generation after mutation"].values, ICdf["# Generation after mutation"].values)),
        'Information': np.concatenate((ICdf["IC"].values, MIdf["MI"].values)),
        "Method": list(np.repeat("Conservation", ICdf.shape[0])) + list(np.repeat("Correlation", MIdf.shape[0])),
        "SEM": np.concatenate((ICdf["SEM_IC"].values,MIdf["SEM_MI"]))
        })

    #Maps Method to color
    colorMap = {"Conservation":"red", "Correlation":"blue"}

    sns.set_style("whitegrid")
    #Creates matlab axis object
    figure = sns.lineplot(
        data=combineddf,x="Generation", y="Information", hue = "Method",  palette = colorMap, ci = None)

    #Adds confidence Interval
    figure.fill_between(ICdf["# Generation after mutation"],ICdf.IC - (ICdf.SEM_IC * 2), ICdf.IC + (ICdf.SEM_IC * 2), alpha = 0.2, color = colorMap["Conservation"])
    figure.fill_between(MIdf["# Generation after mutation"],MIdf.MI - (MIdf.SEM_MI * 2), MIdf.MI + (MIdf.SEM_MI * 2), alpha = 0.2, color = colorMap["Correlation"])
    
    if text == True:
        #Adds labels
        figure.set_title(pltTitle)
        figure.set_xlabel(xLabel)
        figure.set_ylabel(ylabel)
        figure.legend(title = legendLabel)
    else:
        figure.legend([],[],frameon = False)
    

    #Returns figure instance for visualization or download
    return figure.get_figure()


def plotGenericAUCcurve(fileList, text = True, pltTitle = None, xLabel = None, ylabel = None, legendLabel = None):
    """
    Plotting function returning the seaborn figureplot
    Plots MI and IC curves of two different populations. Called by several figure functions
    
    filelist - 2 element list containing strings describing relative path to IC and MI csv files
    text - Boolean indicating whether to include text in the graph
    Returns
    -------
    figureplot

    """
    #Populates dataframe from csv files and creates combined data frame
    ICdf = pd.read_csv(fileList[0])
    MIdf = pd.read_csv(fileList[1])

    #Manipulates combined dataframe for plotting purposes
    combineddf = pd.DataFrame({
        'Generation': np.concatenate((ICdf["Generation"].values, ICdf["Generation"].values)),
        'AUC': np.concatenate((ICdf["AUC"].values, MIdf["AUC"].values)),
        "Method": list(np.repeat("Conservation", ICdf.shape[0])) + list(np.repeat("Correlation", MIdf.shape[0])),
        })

    #Maps Method to color
    colorMap = {"Conservation":"red", "Correlation":"blue"}

    sns.set_style("whitegrid")
    #Creates matlab axis object
    figure = sns.lineplot(
        data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ci = None)

    
    if text == True:
        #Adds labels
        figure.set_title(pltTitle)
        figure.set_xlabel(xLabel)
        figure.set_ylabel(ylabel)
        figure.legend(title = legendLabel)
    else:
        figure.legend([],[],frameon = False)
    

    #Returns figure instance for visualization or download
    return figure.get_figure()

def Figure1(fileList):
    """
    Plotting function returning the seaborn figureplot for Figure 1
    
    filelist - 2 element list containing strings describing relative path to IC and MI csv files
    text - Boolean indicating whether to include text in the graph
    Returns
    -------
    figureplot

    """

    #Legend and Title Information
    pltTitle = "Information loss under mutation for Conservation- and Correlation-encoding populations"
    xLabel = "Generation"
    yLabel = "Information (bits)"
    legendLabel = "Strategy"
    
    #Returns figure created by generic function
    return plotGenericInformationCurve(fileList, text = True, pltTitle = pltTitle, xLabel = xLabel, ylabel = yLabel,legendLabel = legendLabel)




    
def SupplementaryFigure1(fileList):
    """
    

    Parameters
    ----------
    fileList : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    #Legend and Title Information
    pltTitle = "Total Information in various encoding strategies after generations of mutations"
    xLabel = "Generation"
    ylabel = "Information"
    legendLabel = "Strategy                 AUC"

    #Creates combined dataframe from input files
    strategies = ["0.00 IC", "0.25 IC", "0.50 IC", "0.75 IC", "1.00 IC"]
    listDF = [pd.read_csv(x) for x in fileList]
    combineddf = pd.DataFrame({
        'Generation': np.concatenate([x["# Generation after mutation"].values for x in listDF]),
        'Information': np.concatenate([x["TI"].values for x in listDF]),
        'Strategy': np.concatenate([np.repeat(x, listDF[0].shape[0]) for x in strategies]),
        'SEM': np.concatenate([x["SEM_TI"].values for x in listDF])
        })

    #Changes grid style
    sns.set_style("whitegrid")

    #Creates matplotlib axis object containing all TI curves
    colorMap = {"0.00 IC":"darkblue","0.25 IC":"blue", "0.50 IC": "purple", "0.75 IC":"darkred", "1.00 IC":"red"}
    SuppFigure1 = sns.lineplot(
        data=combineddf,x="Generation", y="Information", hue = "Strategy",  palette = colorMap, ci = None)
        
    #Adds confidence intervals for each curve
    for i in range(len(listDF)):
        strategy = listDF[i]
        SuppFigure1.fill_between(strategy["# Generation after mutation"], strategy.TI - (strategy.SEM_TI * 2),strategy.TI + (strategy.SEM_TI * 2), alpha = 0.2, color = colorMap[list(colorMap.keys())[i]])
        
    #Adds labels
    SuppFigure1.set_title(pltTitle)
    SuppFigure1.set_xlabel(xLabel)
    SuppFigure1.set_ylabel(ylabel)

    #Edits labels of legend to include AUC
    AUClist = []
    for i in range(len(listDF)):
        AUClist.append(round(computeAUC(list(listDF[i]["# Generation after mutation"].values), list(listDF[i]["TI"].values), 10),2))

    newlabels = []
    #Alters spacing to align the AUC values in the legend
    for i in  range(len(strategies)):
        strategy = strategies[i]
        
        if len(strategy) == 4:
            space = "       " #normal tab
        else:
            space = "       "; space = space[:(len(space) - (len(strategy) - 4))]
        newlabels.append(strategy + space + str(AUClist[i]))
        
    SuppFigure1.legend(title = legendLabel,labels = newlabels, prop = {'family': 'DejaVu Sans Mono'})

    return SuppFigure1.get_figure()
def SupplementaryFigure2(fileList, subSample = 100):
    """
    

    Parameters
    ----------
    fileList : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    #Legend and Title Information
    pltTitle = "Jensen Shannon Divergence of MI- and IC-encoding populations after mutation"
    xLabel = "Generation"
    yLabel = "Divergence"
    legendLabel = "Strategy             Slope"
    
    #Populates dataframe from csv files and creates combined data frame
    ICdf = pd.read_csv(fileList[0])[0:(subSample + 1)]
    MIdf = pd.read_csv(fileList[1])[0:(subSample + 1)]
    
    ICdf["JSD"] = [(e**x) for x in list(ICdf["JSD"].values)]
    MIdf["JSD"] = [(e**x) for x in list(MIdf["JSD"].values)]
        
    #Manipulates combined dataframe for plotting purposes
    combineddf = pd.DataFrame({
        'Generation': np.concatenate((ICdf["# Generation after mutation"].values, ICdf["# Generation after mutation"].values)),
        'Divergence': np.concatenate((ICdf["JSD"].values, MIdf["JSD"].values)),
        "Method": list(np.repeat("IC", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        #"SEM": np.concatenate((ICdf["SEM_JSD"].values,MIdf["SEM_JSD"]))
        })

    #Maps Method to color
    colorMap = {"IC":"red", "MI":"blue"}

    sns.set_style("whitegrid")
    #Creates matlab axis object
    figure = sns.lineplot(
        data=combineddf,x="Generation", y="Divergence", hue = "Method",  palette = colorMap, ci = None)

    #Adds confidence Interval
   # figure.fill_between(ICdf["# Generation after mutation"],ICdf.JSD - (ICdf.SEM_JSD * 2), ICdf.JSD + (ICdf.SEM_JSD * 2), alpha = 0.2, color = colorMap["IC"])
    #figure.fill_between(MIdf["# Generation after mutation"],MIdf.JSD - (MIdf.SEM_JSD * 2), MIdf.JSD + (MIdf.SEM_JSD * 2), alpha = 0.2, color = colorMap["MI"])
    
    #Adds labels
    figure.set_title(pltTitle)
    figure.set_xlabel(xLabel)
    figure.set_ylabel(yLabel)
    figure.legend(title = legendLabel)
    
    AUClist = []
    
    #For recording AUC
    dfList = [ICdf, MIdf]
    
    for i in range(len(dfList)):
        AUClist.append(round((dfList[i]["JSD"].values[-1] - dfList[i]["JSD"].values[0])/(dfList[i].shape[0]), 4))

    strategies = ["100% IC", "100% MI"]
    newlabels = []
    #Alters spacing to align the AUC values in the legend
    for i in  range(len(dfList)):
        strategy = strategies[i]
        
        if len(strategy) == 4:
            space = "       " #normal tab
        else:
            space = "   ";
        newlabels.append(strategy + space + str(AUClist[i]))
        
    figure.legend(title = legendLabel,labels = newlabels, prop = {'family': 'DejaVu Sans Mono'})
    
    return figure.get_figure()
    
    

    
def Figure2(fileList):
    """
    
    """
    #Creates figure for each pair
    figureList = []

    for populationType in fileList:
        figureList.append(plotGenericInformationCurve(populationType,text = False))
        figure = plt.figure()
        
    #Returns list of figures
    return figureList  

def Figure3(fileList):
    """
    

    Parameters
    ----------
    fileList : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pltTitle = "Information loss under mutation for MI- and IC-encoding populations"
    xLabel = "Generation"
    yLabel = "Information"
    legendLabel = "Strategy"
    
    return plotGenericInformationCurve(fileList,text = True, pltTitle = pltTitle, xLabel = xLabel, ylabel = yLabel, legendLabel = legendLabel)

def Figure4(fileList):
    """
    Plotting function returning the seaborn figureplot for Figure 1
    
    filelist - 2 element list containing strings describing relative path to IC and MI csv files
    text - Boolean indicating whether to include text in the graph
    Returns
    -------
    figureplot

    """

    #Legend and Title Information
    pltTitle = "Area under the precision recall curve after mutation for Conservation- and Correlation-encoding populations"
    xLabel = "Generation"
    yLabel = "AUPRC"
    legendLabel = "Strategy"
    
    #Returns figure created by generic function
    return plotGenericAUCcurve(fileList, text = True, pltTitle = pltTitle, xLabel = xLabel, ylabel = yLabel,legendLabel = legendLabel)