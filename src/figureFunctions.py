# -*- coding: utf-8 -*-
"""
Functions used to generate figures through seaborn. Functions accept file paths
and output figures
"""
import warnings
warnings.filterwarnings("ignore")
import os
import random
from sklearn.metrics import auc
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt    
from statistics import median

def equalizeDFs(listDFs):
    """
    Equalizes the number of rows of inputted dataframes. 
    Sets the number of rows of each dataframe equal to the minimum number 
    of rows in all dataframes in listDF.
    
    listDFs: list
        List of dataframes (IC and MI). Called by figure functions
        
    Returns list of dataframes
    """
    
    #Determines minimum rows
    minimumRows = listDFs[0].shape[0]
    for i in range(len(listDFs)):
        if minimumRows < listDFs[i].shape[0]:
            minimumRows = listDFs[i].shape[0]
    
    return [listDF.iloc[0:minimumRows,] for listDF in listDFs]           
        
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

def centered_subplots(rows,figsize=None):
    
    grid_dim=max(rows)
    grid_shape=(len(rows),2*grid_dim)
    
    if figsize:
        fig = plt.figure(figsize=(figsize))
    else:
        fig = plt.figure(figsize=(2*grid_dim,3*len(rows)))
    
    allaxes=[]
    
    jrow=0
    for row in rows:
        offset=0
        for i in range(row):
            if row<grid_dim:
                offset =grid_dim-row
                
            ax_position=(jrow,2*i+offset)
            ax = plt.subplot2grid(grid_shape, ax_position, fig=fig,colspan=2)
            allaxes.append(ax)
            
        jrow+=1
        
    return fig, allaxes

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
    
    ICdf, MIdf = equalizeDFs([ICdf, MIdf])

    #Manipulates combined dataframe for plotting purposes
    combineddf = pd.DataFrame({
        'Generation': np.concatenate((ICdf["# Generation after mutation"].values, ICdf["# Generation after mutation"].values)),
        'Information': np.concatenate((ICdf["IC"].values, MIdf["MI"].values)),
        "Method": list(np.repeat("IC", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        "SEM": np.concatenate((ICdf["SEM_IC"].values,MIdf["SEM_MI"]))
        })

    #Maps Method to color
    colorMap = {"IC":"red", "MI":"blue"}

    sns.set_style("whitegrid")
    #Creates matlab axis object
    figure = sns.lineplot(
        data=combineddf,x="Generation", y="Information", hue = "Method",  palette = colorMap, ci = None)

    #Adds confidence Interval
    figure.fill_between(ICdf["# Generation after mutation"],ICdf.IC - (ICdf.SEM_IC * 2), ICdf.IC + (ICdf.SEM_IC * 2), alpha = 0.2, color = colorMap["IC"])
    figure.fill_between(MIdf["# Generation after mutation"],MIdf.MI - (MIdf.SEM_MI * 2), MIdf.MI + (MIdf.SEM_MI * 2), alpha = 0.2, color = colorMap["MI"])
    
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
    
    filelist: list
        2 element list with paths to informationDecay csv files for IC/MI populations
        Note: IC path must be specified first and MI path must be specified second

    Returns
    -------
    figureplot

    """

    #Legend and Title Information
    pltTitle = "Information loss under mutation for MI- and IC-encoding populations"
    xLabel = "Generation"
    yLabel = "Information"
    legendLabel = "Strategy"
    
    #Returns figure created by generic function
    return plotGenericInformationCurve(fileList, text = True, pltTitle = pltTitle, xLabel = xLabel, ylabel = yLabel,legendLabel = legendLabel)




    
def SupplementaryFigure1(fileList):
    """
    Plotting function returning the seaborn figureplot for Supplementary Figure 1

    filelist: list
        5 element list with paths to informationDecay csv files for IC/MI populations
        Note: Paths must be specified in the order of 0% IC, 25% IC 50% IC, 75% IC, 
        100% IC

    Returns
    -------
    figureplot

    """
    #Legend and Title Information
    pltTitle = "Total Information in various encoding strategies after generations of mutations"
    xLabel = "Generation"
    ylabel = "Information"
    legendLabel = "Strategy                 AUC"

    #Creates combined dataframe from input files
    strategies = ["0.00 IC", "0.25 IC", "0.50 IC", "0.75 IC", "1.00 IC"]
    listDF = [pd.read_csv(x) for x in fileList]
    
    #Grabs longest 
    listDF = equalizeDFs(listDF)
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

def SupplementaryFigure2(fileList, subSample = None, aucGen = 10):
    """
    Plotting function returning the seaborn figureplot for Supplementary Figure 1

    filelist: list
        2 element list with paths to informationDecay csv files for IC/MI populations
        Note: IC path must be specified first and MI path must be specified second
        
    subsample: int or None
        If specified, subsamples dataframe and cuts the x axis of the plot to generations specified
        
    aucGen: int
        Number of generations to construct linear line representing divergence
        by the first several mutations
        
    Returns
    -------
    figureplot

    """
    #Legend and Title Information
    pltTitle = "Jensen Shannon Divergence of MI- and IC-encoding populations after mutation"
    xLabel = "Mutations"
    yLabel = "Divergence"
    legendLabel = "Strategy             Slope"
    
    
    #Populates dataframe from csv files and creates combined data frame
    ICdf = pd.read_csv(fileList[0])
    MIdf = pd.read_csv(fileList[1])
    
    ICdf, MIdf = equalizeDFs([ICdf,MIdf])
    
    if subSample is not None:
        ICdf = pd.read_csv(fileList[0]).iloc[0:(subSample + 1),]
        MIdf = pd.read_csv(fileList[1]).iloc[0:(subSample + 1),]
    
    #ICdf["# Generation after mutation"] = [log(x + 1, 10) for x in list(ICdf["# Generation after mutation"].values)]
    #MIdf["# Generation after mutation"] = [log(x + 1, 10) for x in list(MIdf["# Generation after mutation"].values)]
    
    #ICdf["JSD"] = [(e**x) for x in list(ICdf["JSD"].values)]
    #MIdf["JSD"] = [(e**x) for x in list(MIdf["JSD"].values)]
        
    #Manipulates combined dataframe for plotting purposes
    combineddf = pd.DataFrame({
        'Generation': np.concatenate((ICdf["# Generation after mutation"].values, ICdf["# Generation after mutation"].values)),
        'Divergence': np.concatenate((ICdf["JSD"].values, MIdf["JSD"].values)),
        "Method": list(np.repeat("IC", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        "SEM": np.concatenate((ICdf["SEM_JSD"].values,MIdf["SEM_JSD"]))
        })
    
    #Maps Method to color
    colorMap = {"IC":"red", "MI":"blue"}
    color2Map = {"IC": "black", "MI":"black"}
    sns.set_style("whitegrid")
    #Creates matlab axis object
    figure = sns.lineplot(
        data=combineddf,x="Generation", y="Divergence", hue = "Method",  palette = colorMap, ci = None)

    #Adds confidence Interval
    figure.fill_between(ICdf["# Generation after mutation"],ICdf.JSD - (ICdf.SEM_JSD * 2), ICdf.JSD + (ICdf.SEM_JSD * 2), alpha = 0.2, color = colorMap["IC"])
    figure.fill_between(MIdf["# Generation after mutation"],MIdf.JSD - (MIdf.SEM_JSD * 2), MIdf.JSD + (MIdf.SEM_JSD * 2), alpha = 0.2, color = colorMap["MI"])
    
    #Adds labels
    figure.set_title(pltTitle)
    figure.set_xlabel(xLabel)
    figure.set_ylabel(yLabel)
    figure.legend(title = legendLabel)
    
    AUClist = []
    
    #For recording AUC
    dfList = [ICdf, MIdf]
    
    for i in range(len(dfList)):
        AUClist.append((dfList[i]["JSD"].values[aucGen] - dfList[i]["JSD"].values[0])/(aucGen))
    
    ICAUCList = [AUClist[0] * x for x in range(len(ICdf["# Generation after mutation"].values))]
    MIAUCList = [AUClist[1] * x for x in range(len(MIdf["# Generation after mutation"].values))]

    combineddf_FirstTenGen = pd.DataFrame({
        'Generation': np.concatenate((ICdf["# Generation after mutation"].values, ICdf["# Generation after mutation"].values)),
        'Divergence': ICAUCList + MIAUCList,
        "Method": list(np.repeat("IC", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        #"SEM": np.concatenate((ICdf["SEM_JSD"].values,MIdf["SEM_JSD"]))
        })
    
    figure = sns.lineplot(
        data=combineddf_FirstTenGen,x="Generation", y="Divergence", hue = "Method",  palette = color2Map, linestyle = '--',ci = None)
    strategies = ["100% IC", "100% MI"]
    newlabels = []
    #Alters spacing to align the AUC values in the legend
    for i in  range(len(dfList)):
        strategy = strategies[i]
        
        if len(strategy) == 4:
            space = "       " #normal tab
        else:
            space = "   ";
        newlabels.append(strategy + space + str(round(AUClist[i],4)))
        
    figure.legend(title = legendLabel,labels = newlabels, prop = {'family': 'DejaVu Sans Mono'})
    
    return figure.get_figure()
    
    

    
def Figure2(fileList):
    """
    Plotting function returning the seaborn figureplot for Figure 2
    
    filelist: list
        4 element list with paths to informationDecay csv files for IC/MI populations
        Note: Inputted files must follow the following order - ["0.5 bits", "1 bit", "1.5 bits", "2 bits"]
            
    Returns
    -------
    figureplot
    
    """
    pltTitle = "Information loss under mutation for MI- and IC-encoding populations"
    xLabel = "Generation"
    yLabel = "AUPRC"
    legendLabel = "Strategy"
    

    #Creates figure for each pair
    figureList = []
    fig,ax = plt.subplots(ncols = 2, nrows = 2) 
    
    Strategies = ["0.5 bits", "1 bit", "1.5 bits", "2 bits"]
    
    for i in range(len(fileList)):
        populationType = fileList[i]
        row = i // 2
        col = i % 2
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
        ICdf = pd.read_csv(populationType[0])
        MIdf = pd.read_csv(populationType[1])
    
        #Manipulates combined dataframe for plotting purposes
        combineddf = pd.DataFrame({
            'Generation': np.concatenate((ICdf["# Generation after mutation"].values, ICdf["# Generation after mutation"].values)),
            'Information': np.concatenate((ICdf["IC"].values, MIdf["MI"].values)),
            "Method": np.concatenate((np.array(np.repeat("IC", ICdf.shape[0])), np.array(np.repeat("MI", MIdf.shape[0])))),
            "SEM": np.concatenate((ICdf["SEM_IC"].values,MIdf["SEM_MI"].values))
            })

    
        #Maps Method to color
        colorMap = {"IC":"red", "MI":"blue"}
    
        sns.set_style("whitegrid")
        #Creates matlab axis object
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="Information", hue = "Method",  palette = colorMap, ax= ax[row,col], errorbar = None)
        
        
        #Adds confidence Interval
        figure.fill_between(ICdf["# Generation after mutation"],ICdf.IC - (ICdf.SEM_IC * 2), ICdf.IC + (ICdf.SEM_IC * 2), alpha = 0.2, color = colorMap["IC"])
        figure.fill_between(MIdf["# Generation after mutation"],MIdf.MI - (MIdf.SEM_MI * 2), MIdf.MI + (MIdf.SEM_MI * 2), alpha = 0.2, color = colorMap["MI"])
        figure.set_title(""); figure.set_xlabel(""); figure.set_ylabel("")
        
        #Position text above plot
        ax[row, col].text(500 - len(Strategies[i]), combineddf["Information"][0] * 1.1, Strategies[i])
        if (i != 3):
            figure.legend([],[],frameon = False)
        else:
            fig.suptitle(pltTitle)
            fig.supxlabel(xLabel)
            fig.supylabel(yLabel)
            handles, labels = ax[1,1].get_legend_handles_labels()
            fig.legend(handles,labels,loc = "upper right")
            figure.legend([],[],frameon = False)
            fig.tight_layout()
            
    #Returns list of figures
    return fig

def Figure3(fileList):
    """
    Plotting function returning the seaborn figureplot for Figure 3
    
    filelist: list
        2 element list with paths to informationDecay csv files for IC/MI populations
        Note: IC path must be specified first and MI path must be specified second

    Returns
    -------
    figureplot

    """
    pltTitle = "Information loss under mutation for MI- and IC-encoding populations"
    xLabel = "Generation"
    yLabel = "Information"
    legendLabel = "Strategy"
    
    return plotGenericInformationCurve(fileList,text = True, pltTitle = pltTitle, xLabel = xLabel, ylabel = yLabel, legendLabel = legendLabel)

def Figure4(fileList):
    """
    Plotting function returning the seaborn figureplot for Figure 4
    
    filelist: list
        4 element list with paths to genomeSearch csv files for IC/MI populations
        Note: Inputted files must follow the following order - ["10 bp", "12 bp", "14 bp", "16 bp"]

    Returns
    -------
    figureplot
    """

    #Legend and Title Information
    pltTitle = "AUPRC after mutation for MI- and IC-encoding populations"
    xLabel = "Generation"
    yLabel = "AUPRC"
    legendLabel = "Strategy"
    
    dfList = []
    for collection in fileList:
        for file in collection:
            dfList.append(pd.read_csv(file))
            
    dfList = equalizeDFs(dfList)
    #Creates figure for each pair
    figureList = []
    fig,ax = plt.subplots(ncols = 2, nrows = 2) 
    
    Strategies = ["10 bp", "12 bp", "14 bp", "16 bp"]
    for i in range(len(fileList)):
        populationType = fileList[i]
        dfNum = i * 2
        row = i // 2
        col = i % 2
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
        ICdf = pd.read_csv(populationType[0])
        MIdf = pd.read_csv(populationType[1])
        

        #Manipulates combined dataframe for plotting purposes
        combineddf = pd.DataFrame({
            'Generation': np.concatenate((ICdf["Generation"].values, MIdf["Generation"].values)),
            'AUC': np.concatenate((ICdf["AUC"].values, MIdf["AUC"].values)),
            "Method": list(np.repeat("IC", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
            #"SEM": np.concatenate((ICdf["SEM_IC"].values,MIdf["SEM_MI"]))
            })
    
        #Maps Method to color
        colorMap = {"IC":"red", "MI":"blue"}
    
        #sns.set_style("whitegrid")
        #Creates matlab axis object
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ax= ax[row,col], errorbar = None)
        
        ax[row,col].yaxis.set_major_locator(plt.MaxNLocator(4))
        #Adds confidence Interval
        #figure.fill_between(ICdf["Generation"],ICdf.IC - (ICdf.SEM_IC * 2), ICdf.IC + (ICdf.SEM_IC * 2), alpha = 0.2, color = colorMap["IC"])
        #figure.fill_between(MIdf["Generation"],MIdf.MI - (MIdf.SEM_MI * 2), MIdf.MI + (MIdf.SEM_MI * 2), alpha = 0.2, color = colorMap["MI"])
        figure.set_title(""); figure.set_xlabel(""); figure.set_ylabel("")
        
        ax[row, col].text(50, combineddf["AUC"][0] * 1.1, Strategies[i])
        if (i != 3):
            figure.legend([],[],frameon = False)
        else:
            fig.suptitle(pltTitle)
            fig.supxlabel(xLabel)
            fig.supylabel(yLabel)
            handles, labels = ax[1,1].get_legend_handles_labels()
            fig.legend(handles,labels,loc = "upper right")
            figure.legend([],[],frameon = False)
            fig.tight_layout()
            
    #Returns list of figures
    return fig

def Figure5(fileList, fullYaxis = False):
    """
    Plotting function returning the seaborn figureplot for Figure 5
    
    filelist: list
        3 element list with paths to genomeSearch csv files for IC/MI populations
        Note: Inputted files must follow the following order - ["10 bp", "12 bp", "14 bp", "16 bp"]

    Returns
    -------
    figureplot

    """

    #Legend and Title Information
    pltTitle = "AUPRC after mutation for MI- and IC-encoding populations"
    xLabel = "Generation"
    yLabel = "AUPRC"
    legendLabel = "Strategy"
    
    dfList = []
    for collection in fileList:
        for file in collection:
            dfList.append(pd.read_csv(file))
            
    dfList = equalizeDFs(dfList)
    #Creates figure for each pair
    figureList = []
    #fig,ax = plt.subplots(ncols = 2, nrows = 2) 
    fig, allaxes = centered_subplots([2,1], [6,4])
    
    Strategies = ["TI: 4.8 bit", "TI: 6 bit", "TI: 7.2 bit"]
    for i in range(len(fileList)):
        populationType = fileList[i]
        dfNum = i * 2
        row = i // 2
        col = i % 2
        print(row,col)
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
        ICdf = pd.read_csv(populationType[0])
        MIdf = pd.read_csv(populationType[1])
        

        #Manipulates combined dataframe for plotting purposes
        combineddf = pd.DataFrame({
            'Generation': np.concatenate((ICdf["Generation"].values, MIdf["Generation"].values)),
            'AUC': np.concatenate((ICdf["AUC"].values, MIdf["AUC"].values)),
            "Method": list(np.repeat("IC", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
            #"SEM": np.concatenate((ICdf["SEM_IC"].values,MIdf["SEM_MI"]))
            })
    
        #Maps Method to color
        colorMap = {"IC":"red", "MI":"blue"}
    
        #sns.set_style("whitegrid")
        #Creates matlab axis object
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ax= allaxes[i], errorbar = None)
        
        
        #Adds confidence Interval
        #figure.fill_between(ICdf["Generation"],ICdf.AUC - (ICdf.SEM_IC * 2), ICdf.IC + (ICdf.SEM_IC * 2), alpha = 0.2, color = colorMap["IC"])
        #figure.fill_between(MIdf["Generation"],MIdf.MI - (MIdf.SEM_MI * 2), MIdf.MI + (MIdf.SEM_MI * 2), alpha = 0.2, color = colorMap["MI"])
        figure.set_title(""); figure.set_xlabel(""); figure.set_ylabel(""); 
        
        if fullYaxis:
            figure.set_ylim([0,1])
            allaxes[i].text(0, 1*1.05, Strategies[i])
        else:
            figure.set_ylim([0,combineddf["AUC"][0] * 1.1])
            allaxes[i].text(0, combineddf["AUC"][0] * 1.15, Strategies[i])  

        allaxes[i].yaxis.set_major_locator(plt.MaxNLocator(4))

        if (i <= 1):
            figure.legend([],[],frameon = False)
        else:
            fig.suptitle(pltTitle)
            fig.supxlabel(xLabel)
            fig.supylabel(yLabel)
            handles, labels = allaxes[i].get_legend_handles_labels()
            fig.legend(handles,labels,loc = "upper right")
            figure.legend([],[],frameon = False)
            #fig.tight_layout()
            plt.subplots_adjust(hspace = 0.3, wspace = 0.4)
    
    
    """Make sure to scale y axis starting from 0"""
    #Returns list of figures
    return fig

def Figure6(subdir):
    """
    Plotting function returning the seaborn figureplot for Figure 6
    
    subdir: str
        Path to folder containing folders of competition experiment
        
    Returns
    -------
    figureplot

    """

    #Grabs file paths for each competition experiment
    competitionPaths = [x[0].replace('\\','/') + "/runCompetition.csv" for x in os.walk(subdir)]
    competitionPaths.pop(0) #Removes nonsubdirectory
    competitionDFs = [pd.read_csv(x) for x in competitionPaths]
    
    #Adds key to every run
    for i in range(len(competitionDFs)):
        competitionDFs[i]["Run Tag"] = [f"Run {i}" for k in range(competitionDFs[i].shape[0])]
    
    #Constructs average DF for average curve
    averageDF = competitionDFs[0]
    averageDF["Run Tag"] = ["Average" for k in range(competitionDFs[0].shape[0])]
    for i in range(competitionDFs[0].shape[0]):
        averageDF["MI_80 org frequency"][i] = median([list(df["MI_80 org frequency"])[i] for df in competitionDFs])
    
    #Constructs cumulative pandas DF with all data
    cumulativeDF = pd.concat(competitionDFs)
    
    ax = sns.lineplot(data = cumulativeDF, x = "Generation", y = "MI_80 org frequency", 
                      units = "Run Tag", estimator = None, alpha = 0.5, linewidth = 0.5,
                      color = "grey")
    ax = sns.lineplot(data = averageDF, x = "Generation", y = "MI_80 org frequency",
                      color = "blue")
    return ax.get_figure()

def Figure7(fileList, fullYaxis = False):
    """
    Plotting function returning the seaborn figureplot for Figure 5
    
    filelist: list
        5 element list with paths to genomeSearch csv files for IC/MI populations
        Note: Inputted files must follow the following order - ["8","10 bp", "12 bp", "15 bp", "20 bp"]

    Returns
    -------
    figureplot

    """

    #Legend and Title Information
    pltTitle = "AUPRC after mutation for MI- and IC-encoding populations"
    xLabel = "Generation"
    yLabel = "AUPRC"
    legendLabel = "Strategy"
    
    dfList = []
    for collection in fileList:
        for file in collection:
            dfList.append(pd.read_csv(file))
            
    dfList = equalizeDFs(dfList)
    #Creates figure for each pair
    figureList = []
    #fig,ax = plt.subplots(ncols = 2, nrows = 2) 
    fig, allaxes = centered_subplots(centered_subplots([3,2], [6,4]))
    
    Strategies = ["8bp","10bp", "12bp", "14bp","16bp"]
    for i in range(len(fileList)):
        populationType = fileList[i]
        dfNum = i * 2
        row = i // 2
        col = i % 2
        print(row,col)
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
        ICdf = pd.read_csv(populationType[0])
        MIdf = pd.read_csv(populationType[1])
        

        #Manipulates combined dataframe for plotting purposes
        combineddf = pd.DataFrame({
            'Generation': np.concatenate((ICdf["Generation"].values, MIdf["Generation"].values)),
            'AUC': np.concatenate((ICdf["AUC"].values, MIdf["AUC"].values)),
            'SEM_AUC': np.concatenate((ICdf["SEM_AUC"].values, MIdf["SEM_AUC"].values)),
            "Method": list(np.repeat("IC", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
            #"SEM": np.concatenate((ICdf["SEM_IC"].values,MIdf["SEM_MI"]))
            })
    
        #Maps Method to color
        colorMap = {"IC":"red", "MI":"blue"}
    
        #sns.set_style("whitegrid")
        #Creates matlab axis object
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ax= allaxes[i], errorbar = None)
        
        
        #Adds confidence Interval
        figure.fill_between(ICdf["Generation"],ICdf.AUC - (ICdf.SEM_AUC * 2), ICdf.AUC + (ICdf.SEM_AUC * 2), alpha = 0.2, color = colorMap["IC"])
        figure.fill_between(MIdf["Generation"],MIdf.AUC - (MIdf.SEM_AUC * 2), MIdf.AUC + (MIdf.SEM_AUC * 2), alpha = 0.2, color = colorMap["MI"])
        figure.set_title(""); figure.set_xlabel(""); figure.set_ylabel(""); 
        
        if fullYaxis:
            figure.set_ylim([0,1])
            allaxes[i].text(0, 1*1.05, Strategies[i])
        else:
            figure.set_ylim([0,combineddf["AUC"][0] * 1.1])
            allaxes[i].text(0, combineddf["AUC"][0] * 1.15, Strategies[i])  

        allaxes[i].yaxis.set_major_locator(plt.MaxNLocator(4))

        if (i <= 1):
            figure.legend([],[],frameon = False)
        else:
            fig.suptitle(pltTitle)
            fig.supxlabel(xLabel)
            fig.supylabel(yLabel)
            handles, labels = allaxes[i].get_legend_handles_labels()
            fig.legend(handles,labels,loc = "upper right")
            figure.legend([],[],frameon = False)
            #fig.tight_layout()
            plt.subplots_adjust(hspace = 0.3, wspace = 0.4)
    
    
    """Make sure to scale y axis starting from 0"""
    #Returns list of figures
    return fig
