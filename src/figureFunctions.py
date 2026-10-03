# -*- coding: utf-8 -*-
"""
Functions used to generate figures through seaborn. Functions accept file paths
and output figures. Each function plots a single panel; see mainFigures.py for the
paper-numbered plates that assemble these panels.
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
from statistics import median, mean
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import copy
import math
import importlib
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from src.Objects.mi_vs_si import MI_VS_SI
import matplotlib.patches as patches

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

def plotGenericInformationCurve(fileList, text=True, pltTitle=None, xLabel=None, ylabel=None,
                                legendLabel=None, axis=None, yticks=5):
    """
    Plotting function returning the seaborn figureplot
    """

    ICdf = pd.read_csv(fileList[0])
    MIdf = pd.read_csv(fileList[1])

    ICdf, MIdf = equalizeDFs([ICdf, MIdf])

    combineddf = pd.DataFrame({
        'Generation': np.concatenate((ICdf["# Generation after mutation"].values,
                                       MIdf["# Generation after mutation"].values)),
        'Information': np.concatenate((ICdf["IC"].values, MIdf["MI"].values)),
        'Method': list(np.repeat("SI", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        'SEM': np.concatenate((ICdf["SEM_IC"].values, MIdf["SEM_MI"].values))
    })

    colorMap = {"SI": "red", "MI": "blue"}

    sns.set_style("whitegrid")
    figure = sns.lineplot(
        data=combineddf, x="Generation", y="Information", hue="Method",
        palette=colorMap, ci=None, ax=axis
    )

    figure.fill_between(
        ICdf["# Generation after mutation"],
        ICdf.IC - (ICdf.SEM_IC * 2),
        ICdf.IC + (ICdf.SEM_IC * 2),
        alpha=0.2,
        color=colorMap["SI"]
    )
    figure.fill_between(
        MIdf["# Generation after mutation"],
        MIdf.MI - (MIdf.SEM_MI * 2),
        MIdf.MI + (MIdf.SEM_MI * 2),
        alpha=0.2,
        color=colorMap["MI"]
    )

    if text:
        figure.set_title(pltTitle)
        figure.set_xlabel(xLabel)
        figure.set_ylabel(ylabel)
        figure.legend(title=legendLabel)
    else:
        figure.legend([], [], frameon=False)

    if yticks is not None:
        figure.yaxis.set_major_locator(MaxNLocator(nbins=yticks))

    return figure.get_figure()

def informationDecayFig(fileList, axis=None):
    """
    Plotting function returning the seaborn figureplot of information decay
    over mutations for a PI- and an MI-encoding population

    filelist: list
        2 element list with paths to informationDecay csv files for PI/MI populations
        Note: IC path must be specified first and MI path must be specified second

    Axis: axes matploblib object or None
        Used for main figures
    Returns
    -------
    figureplot

    """

    pltTitle = ""
    xLabel = "Mutations"
    yLabel = "Information (bits)"
    legendLabel = "Strategy"

    if axis:
        fig = plotGenericInformationCurve(fileList, text=True, pltTitle=pltTitle, xLabel=xLabel, ylabel=yLabel, legendLabel=legendLabel, axis=axis)
        ax = axis
    else:
        fig = plotGenericInformationCurve(fileList, text=True, pltTitle=pltTitle, xLabel=xLabel, ylabel=yLabel, legendLabel=legendLabel)
        ax = fig.axes[0]

    ax.set_xlim(0, 500)

    # Remove legend (shared legend handled in MainFigure2)
    legend = ax.get_legend()
    if legend:
        legend.remove()

    # Read data for inset slope computation
    ICdf = pd.read_csv(fileList[0])
    MIdf = pd.read_csv(fileList[1])
    gen_col = "# Generation after mutation"

    ic_slope = (ICdf["IC"].iloc[1] - ICdf["IC"].iloc[0]) / (ICdf[gen_col].iloc[1] - ICdf[gen_col].iloc[0])
    mi_slope = (MIdf["MI"].iloc[1] - MIdf["MI"].iloc[0]) / (MIdf[gen_col].iloc[1] - MIdf[gen_col].iloc[0])

    ic_x_end = min(-ICdf["IC"].iloc[0] / ic_slope, 100)
    mi_x_end = min(-MIdf["MI"].iloc[0] / mi_slope, 100)

    # Create inset zoomed into first 100 generations with slope lines
    axins = inset_axes(ax, width="50%", height="50%", loc='upper right')

    ICdf_sub = ICdf[ICdf[gen_col] <= 100]
    MIdf_sub = MIdf[MIdf[gen_col] <= 100]

    colorMap = {"SI": "red", "MI": "blue"}
    axins.plot(ICdf_sub[gen_col], ICdf_sub["IC"], color=colorMap["SI"])
    axins.plot(MIdf_sub[gen_col], MIdf_sub["MI"], color=colorMap["MI"])

    x_ic = np.array([0, ic_x_end])
    x_mi = np.array([0, mi_x_end])
    axins.plot(x_ic, ICdf["IC"].iloc[0] + ic_slope * x_ic, color="red", linestyle=":", linewidth=1.5)
    axins.plot(x_mi, MIdf["MI"].iloc[0] + mi_slope * x_mi, color="blue", linestyle=":", linewidth=1.5)

    mi_val_at_100 = MIdf.loc[MIdf[gen_col] == 100, "MI"].values[0]
    axins.set_xlim(0, 100)
    axins.set_ylim(mi_val_at_100, 2)
    axins.set_xticks([0, 25, 50, 75, 100])
    axins.set_yticks([t for t in [0, 0.5, 1, 1.5, 2] if t >= mi_val_at_100])
    axins.set_ylabel("")
    axins.set_xlabel("")

    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    return fig


def informationDistributionFig(popFile):
    """
    Plotting function returning the figureplot for Supplementary Figure S1

    popFile: str
        Path to a motifPopulation json file. Assumes 2 column motifs, since the
        PI/IC of the first column is plotted against the PI/IC of the second

    Returns
    -------
    figureplot

    """
    import json

    with open(popFile, 'r') as f:
        pop_data = json.load(f)

    population = pop_data["populationList"]
    goal_TI = pop_data["GoalIC"] + pop_data["GoalMI"]

    col1_IC = []
    col2_IC = []

    for organism_bs in population:
        N = len(organism_bs)
        L = len(organism_bs[0])
        for col in range(L):
            counts = {"G": 0, "C": 0, "A": 0, "T": 0}
            for bs in organism_bs:
                counts[bs[col]] += 1
            Hx = 0
            for base, count in counts.items():
                freq = count / N
                if freq > 0:
                    Hx -= freq * np.log2(freq)
            ic_col = 2 - Hx
            if col == 0:
                col1_IC.append(ic_col)
            else:
                col2_IC.append(ic_col)

    fig = plt.figure(figsize=(6, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                          hspace=0, wspace=0)

    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    ax_main.scatter(col1_IC, col2_IC, color='red', alpha=0.6, edgecolors='darkred', linewidths=0.5)
    ax_main.plot([0, goal_TI], [goal_TI, 0], 'k--', linewidth=1)
    ax_main.set_xlim(0, 2)
    ax_main.set_ylim(0, 2)
    ax_main.set_xlabel("SI Position 1 (bits)")
    ax_main.set_ylabel("SI Position 2 (bits)")
    ax_main.set_aspect('equal')

    ax_top.hist(col1_IC, bins=10, color='grey', edgecolor='black', linewidth=0.5)
    ax_top.axis('off')

    ax_right.hist(col2_IC, bins=10, color='grey', edgecolor='black', linewidth=0.5,
                  orientation='horizontal')
    ax_right.axis('off')

    fig.add_subplot(gs[0, 1]).axis('off')

    return fig


def mixedPopulationDecayFig(fileList):
    """
    Plotting function returning the seaborn figureplot of information decay
    for populations mixing the PI and MI encoding strategies

    filelist: list
        5 element list with paths to informationDecay csv files for PI/MI populations
        Note: Paths must be specified in the order of 0% PI, 25% PI, 50% PI, 75% PI, 
        100% IC

    Returns
    -------
    figureplot

    """
    #Legend and Title Information
    pltTitle = ""
    xLabel = "Mutations"
    ylabel = "Information (bits)"
    legendLabel = "Strategy                 AUC"

    #Creates combined dataframe from input files
    strategies = ["0.00 SI", "0.25 SI", "0.50 SI", "0.75 SI", "1.00 SI"]
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
    colorMap = {"0.00 SI":"blue","0.25 SI":"#7F40BF", "0.50 SI": "purple", "0.75 SI":"#BF4040", "1.00 SI":"red"}
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

def divergenceFig(fileList, subSample = None, aucGen = 50, axis=None, ylim=1.4, inset_xmax = 200, inset_ymax = 0.3):
    """
    Plotting function returning the seaborn figureplot of Jensen-Shannon divergence
    over mutations for a PI- and an MI-encoding population

    filelist: list
        2 element list with paths to informationDecay csv files for PI/MI populations
        Note: IC path must be specified first and MI path must be specified second
        
    subsample: int or None
        If specified, subsamples dataframe and cuts the x axis of the plot to generations specified
        
    aucGen: int
        Number of generations to construct linear line representing divergence
        by the first several mutations
    
    Axis: axes matploblib object or None
        Used for main figures
    
    ylim: int
        Maximum value of main y axis
        
    inset_xmax: int
        Inset zooms into from x = 0 to x = 250 
        
    inset_ymax: int
        Inset zooms into from y = 0 to Y = 250 
        
    
    
    Returns
    -------
    figureplot

    """
    #Legend and Title Information
    pltTitle = ""
    xLabel = "Mutations"
    yLabel = "Divergence"
    legendLabel = "Strategy"
    
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
        "Method": list(np.repeat("SI", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        "SEM": np.concatenate((ICdf["SEM_JSD"].values,MIdf["SEM_JSD"]))
        })

    #Maps Method to color
    colorMap = {"SI":"red", "MI":"blue"}
    color2Map = {"SI": "black", "MI":"black"}
    sns.set_style("whitegrid")

    if axis:
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="Divergence", hue = "Method",  palette = colorMap, ci = None,ax=axis)
        axis.set_ylim(0,ylim)
    else:
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="Divergence", hue = "Method",  palette = colorMap, ci = None)

    #Adds confidence Interval
    figure.fill_between(ICdf["# Generation after mutation"],ICdf.JSD - (ICdf.SEM_JSD * 2), ICdf.JSD + (ICdf.SEM_JSD * 2), alpha = 0.2, color = colorMap["SI"])
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
        "Method": list(np.repeat("SI", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        })


    if axis:
        # Create inset
        axins = inset_axes(axis, width="60%", height="60%", loc='upper right')
        x1, x2 = 0, inset_xmax
        y1, y2 = 0, inset_ymax

        #Plots
        sns.lineplot(
            data=combineddf,x="Generation", y="Divergence", hue = "Method",  palette = colorMap, ci = None,ax=axins, legend=False)

        sns.lineplot(
            data=combineddf_FirstTenGen,
            x="Generation",
            y="Divergence",
            hue="Method",
            style="Method",
            dashes={"SI": (5, 3), "MI": (1, 1)},
            palette=color2Map,
            ci=None,
            ax=axins,
            legend=True
        )

        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)

        mark_inset(axis, axins, loc1=2, loc2=4, fc="none", ec="0.5")

        axins.set_ylabel("")
        axins.set_xlabel("")

        # Get updated handles for legend
        handles, labels = axins.get_legend_handles_labels()
        axins.legend(handles=handles, labels=[
            f"SI {round(AUClist[0], 4)}",
            f"MI {round(AUClist[1], 4)}"],
            title="Divergence Rate\n(bits/mutation)", title_fontsize="medium",
    labelspacing=1.0,
    handlelength=4.0)

    else:
        figure = sns.lineplot(
            data=combineddf_FirstTenGen,x="Generation", y="Divergence", hue = "Method",  palette = color2Map, linestyle = '--',ci = None)

    strategies = ["100% SI", "100% MI"]
    newlabels = []
    #Alters spacing to align the AUC values in the legend
    for i in  range(len(dfList)):
        strategy = strategies[i]
        
        if len(strategy) == 4:
            space = "       " #normal tab
        else:
            space = "       ";
        #newlabels.append(strategy + space + str(round(AUClist[i],4)))
        newlabels.append(strategy)
        
    #figure.legend(title = legendLabel,labels = newlabels, title_fontsize=10, prop = {'family': 'DejaVu Sans Mono', 'size':9}, 
                 # columnspacing=1.5, handlelength=1, loc="upper left", borderaxespad=0.1,
                 # bbox_to_anchor=(0, 0.75,1, 0.2))
    
    figure.legend(title = legendLabel,labels = newlabels, prop = {'family': 'DejaVu Sans Mono'}, loc="upper left")
    
    return figure.get_figure()
    
    

    
def informationDecayByTIFig(fileList):
    """
    Plotting function returning the multipane seaborn figureplot of information decay
    across total information (TI) values

    filelist: list
        4 element list with paths to informationDecay csv files for PI/MI populations
        Note: Inputted files must follow the following order - ["0.5 bits", "1 bit", "1.5 bits", "2 bits"]

    Returns
    -------
    figureplot

    """
    pltTitle = ""
    xLabel = "Mutations"
    yLabel = "Information (bits)"
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
        Plots MI and IC/PI curves of two different populations. Called by several figure functions
        
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
            "Method": np.concatenate((np.array(np.repeat("SI", ICdf.shape[0])), np.array(np.repeat("MI", MIdf.shape[0])))),
            "SEM": np.concatenate((ICdf["SEM_IC"].values,MIdf["SEM_MI"].values))
            })


        #Maps Method to color
        colorMap = {"SI":"red", "MI":"blue"}

        sns.set_style("whitegrid")
        #Creates matlab axis object
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="Information", hue = "Method",  palette = colorMap, ax= ax[row,col], errorbar = None)


        #Adds confidence Interval
        figure.fill_between(ICdf["# Generation after mutation"],ICdf.IC - (ICdf.SEM_IC * 2), ICdf.IC + (ICdf.SEM_IC * 2), alpha = 0.2, color = colorMap["SI"])
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

def residualInformationFig(fileList, axis=None, L=12, N=128):
    """
    Plotting function returning the seaborn figureplot of information decay.
    Adds dashed and dotted horizontal lines for expected residual PI and MI.

    Parameters
    ----------
    fileList : list[str]
        Two-element list with paths to informationDecay CSV files for PI/MI populations.
        (IC/PI path must be first; MI path second)
    axis : matplotlib.axes.Axes, optional
        Axis to draw on (used for main figures). If None, a new figure is created.
    L : int
        Length of the motif (number of positions in the binding site).
    N : int
        Number of binding sites (sample size used to estimate information).

    Returns
    -------
    matplotlib.figure.Figure or matplotlib.axes.Axes
        The plot or the axis, depending on how the function is called.
    """
    pltTitle = ""
    xLabel = "Mutations"
    yLabel = "Information (bits)"
    legendLabel = "Strategy"

    # Compute expected residual values for DNA (S=4)
    S = 4
    expected_PI = (L * (S - 1)) / (2 * N * math.log(2))          # bits
    expected_MI = ((L - 1) * (S - 1)**2) / (2 * N * math.log(2)) # bits

    # Create or use the axis
    if axis is not None:
        fig = plotGenericInformationCurve(
            fileList, text=True, pltTitle=pltTitle,
            xLabel=xLabel, ylabel=yLabel,
            legendLabel=legendLabel, axis=axis
        )
    else:
        fig = plotGenericInformationCurve(
            fileList, text=True, pltTitle=pltTitle,
            xLabel=xLabel, ylabel=yLabel,
            legendLabel=legendLabel
        )
    ax = fig.axes[0]
    # Add horizontal lines for expected PI and MI
    ax.axhline(expected_PI, color="red", linestyle="--",
               linewidth=1.2, label=f"_nolegend_")
    ax.axhline(expected_MI, color="blue", linestyle=":",
               linewidth=1.2, label=f"_nolegend_")

    # Optional: ensure legend includes both lines
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, frameon=False)

    return fig
def genomeSearchByLengthFig(fileList):
    """
    Plotting function returning the multipane seaborn figureplot of genome search
    performance (AUPRC) across motif lengths

    filelist: list
        4 element list with paths to genomeSearch csv files for PI/MI populations
        Note: Inputted files must follow the following order - ["10 bp", "12 bp", "14 bp", "16 bp"]

    Returns
    -------
    figureplot
    """

    #Legend and Title Information
    pltTitle = ""
    xLabel = "Mutations"
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

        #Populates dataframe from csv files and creates combined data frame
        ICdf = pd.read_csv(populationType[0])
        MIdf = pd.read_csv(populationType[1])


        #Manipulates combined dataframe for plotting purposes
        combineddf = pd.DataFrame({
            'Generation': np.concatenate((ICdf["Generation"].values, MIdf["Generation"].values)),
            'AUC': np.concatenate((ICdf["AUC"].values, MIdf["AUC"].values)),
            "Method": list(np.repeat("SI", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
            })

        #Maps Method to color
        colorMap = {"SI":"red", "MI":"blue"}

        #Creates matlab axis object
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ax= ax[row,col], errorbar = None)

        ax[row,col].yaxis.set_major_locator(plt.MaxNLocator(4))
        figure.set_title(""); figure.set_xlabel(""); figure.set_ylabel("")

        ax[row, col].text(50, combineddf["AUC"][0] * 1.1, Strategies[i])
        if (i != 3):
            figure.legend([],[],frameon = False)
        else:
            fig.supxlabel(xLabel)
            fig.supylabel(yLabel)
            handles, labels = ax[1,1].get_legend_handles_labels()
            fig.legend(handles,labels,loc = "upper right")
            figure.legend([],[],frameon = False)
            fig.tight_layout()
            
    #Returns list of figures
    return fig

def genomeSearchFig(fileList, axis=None):
    """
    Plotting function returning the seaborn figureplot of genome search
    performance (AUPRC) over mutations

    filelist: list
        4 element list with paths to genomeSearch csv files for PI/MI populations
        Note: Inputted files must follow the following order - ["10 bp", "12 bp", "14 bp", "16 bp"]
        
    Axis: axes matploblib object or None
        Used for main figures
    Returns
    -------
    figureplot
    """

    #Legend and Title Information
    #pltTitle = "AUPRC after mutation for MI- and PI-encoding populations"
    pltTitle=""
    xLabel = "Mutations"
    yLabel = "AUPRC"
    legendLabel = "Strategy"
    
    #Populates dataframe from csv files and creates combined data frame
    ICdf = pd.read_csv(fileList[0])
    MIdf = pd.read_csv(fileList[1])
    
    ICdf, MIdf = equalizeDFs([ICdf, MIdf])

    #Manipulates combined dataframe for plotting purposes
    combineddf = pd.DataFrame({
        'Generation': np.concatenate((ICdf["Generation"].values, MIdf["Generation"].values)),
        'AUC': np.concatenate((ICdf["AUC"].values, MIdf["AUC"].values)),
        "Method": list(np.repeat("SI", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
        "SEM_AUC": np.concatenate((ICdf["SEM_AUC"].values,MIdf["SEM_AUC"]))
        })

    #Maps Method to color
    colorMap = {"SI":"red", "MI":"blue"}

    sns.set_style("whitegrid")

    if axis:
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ci = None, ax=axis)
    else:
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ci = None)
    #Adds confidence Interval
    figure.fill_between(ICdf["Generation"],ICdf.AUC - (ICdf.SEM_AUC * 2), ICdf.AUC + (ICdf.SEM_AUC * 2), alpha = 0.2, color = colorMap["SI"])
    figure.fill_between(MIdf["Generation"],MIdf.AUC - (MIdf.SEM_AUC * 2), MIdf.AUC + (MIdf.SEM_AUC * 2), alpha = 0.2, color = colorMap["MI"])
    

    #Adds labels
    figure.set_title(pltTitle)
    figure.set_xlabel(xLabel)
    figure.set_ylabel(yLabel)
    figure.legend(title = legendLabel)

    

    #Returns figure instance for visualization or download
    return figure.get_figure()

def genomeSearchByTIFig(fileList, fullYaxis = False):
    """
    Plotting function returning the multipane seaborn figureplot of genome search
    performance (AUPRC) across total information (TI) values

    filelist: list
        3 element list with paths to genomeSearch csv files for IC/MI populations
        Note: Inputted files must follow the following order - ["10 bp", "12 bp", "14 bp", "16 bp"]

    Returns
    -------
    figureplot

    """

    #Legend and Title Information
    pltTitle = ""
    xLabel = "Mutations"
    yLabel = "AUPRC"
    legendLabel = "Strategy"

    dfList = []
    for collection in fileList:
        for file in collection:
            dfList.append(pd.read_csv(file))

    dfList = equalizeDFs(dfList)
    #Creates figure for each pair
    figureList = []
    fig, allaxes = centered_subplots([2,1], [6,4])

    Strategies = ["TI: 4.8 bit", "TI: 6 bit", "TI: 7.2 bit"]
    for i in range(len(fileList)):
        populationType = fileList[i]
        dfNum = i * 2
        row = i // 2
        col = i % 2

        #Populates dataframe from csv files and creates combined data frame
        ICdf = pd.read_csv(populationType[0])
        MIdf = pd.read_csv(populationType[1])


        #Manipulates combined dataframe for plotting purposes
        combineddf = pd.DataFrame({
            'Generation': np.concatenate((ICdf["Generation"].values, MIdf["Generation"].values)),
            'AUC': np.concatenate((ICdf["AUC"].values, MIdf["AUC"].values)),
            "Method": list(np.repeat("SI", ICdf.shape[0])) + list(np.repeat("MI", MIdf.shape[0])),
            })

        #Maps Method to color
        colorMap = {"SI":"red", "MI":"blue"}

        #Creates matlab axis object
        figure = sns.lineplot(
            data=combineddf,x="Generation", y="AUC", hue = "Method",  palette = colorMap, ax= allaxes[i], errorbar = None)

        figure.set_title(""); figure.set_xlabel(""); figure.set_ylabel("")

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
            fig.supxlabel(xLabel)
            fig.supylabel(yLabel)
            handles, labels = allaxes[i].get_legend_handles_labels()
            fig.legend(handles,labels,loc = "upper right")
            figure.legend([],[],frameon = False)
            plt.subplots_adjust(hspace = 0.3, wspace = 0.4)
    
    
    """Make sure to scale y axis starting from 0"""
    #Returns list of figures
    return fig


def competitionFig(subdir = "./Output/resultsCompetition/", subfolder="competitionpair",fileName = "runCompetition.csv", target_len=1000, axis=None):
    """
    Plotting function returning the seaborn figureplot for competition runs
    
    subdir: str
        Path to folder containing folders of competition experiment
        
    Axis: axes matploblib object or None
        Used for main figures
    Returns
    -------
    figureplot

    """

    #Grabs file paths for each competition experiment
    competitionDFs = [pd.read_csv(f"{subdir}/{subfolder}_{i}/{fileName}") for i in range(100)]
    
    # Step 3: Pad each DataFrame with its last row
    extended_dfs = []
    for df in competitionDFs:
        current_len = df.shape[0]
        if current_len < target_len:
            last_row = df.iloc[[-1]].copy()  # select last row as a DataFrame
            padding = pd.concat([last_row] * (target_len - current_len), ignore_index=True)
            extended_df = pd.concat([df, padding], ignore_index=True)
        else:
            extended_df = df.reset_index(drop=True)
        extended_dfs.append(extended_df)
    competitionDFs = extended_dfs
    print(competitionDFs[0])
    for i in range(len(competitionDFs)):
        competitionDFs[i]["Run Tag"] = [f"Run {i}" for k in range(competitionDFs[i].shape[0])]
        competitionDFs[i]["MI_70 org frequency"] = [1 - freq for freq in competitionDFs[i]["MI_70 org frequency"]]
    
    #Constructs average DF for average curve
    averageDF = copy.deepcopy(competitionDFs[0])
    averageDF["Generation"] = range(target_len)
    averageDF["Run Tag"] = ["Average" for k in range(competitionDFs[0].shape[0])]
    for i in range(target_len):
        averageDF["MI_70 org frequency"][i] = median([list(df["MI_70 org frequency"])[i] for df in competitionDFs])
        
    #Constructs cumulative pandas DF with all data
    cumulativeDF = pd.concat(competitionDFs)
    
    if axis:
        ax = sns.lineplot(data = cumulativeDF, x = "Generation", y = "MI_70 org frequency", 
                          units = "Run Tag", estimator = None, alpha = 0.5, linewidth = 0.5,
                          color = "grey",ax=axis)
        ax = sns.lineplot(data = averageDF, x = "Generation", y = "MI_70 org frequency",
                          color = "blue",ax=axis)
        ax.set_ylabel("SI strategists' frequency")
        #ax.set_xlim(averageDF["Generation"].min(), averageDF["Generation"].max())
    else:
        ax = sns.lineplot(data = cumulativeDF, x = "Generation", y = "MI_70 org frequency", 
                          units = "Run Tag", estimator = None, alpha = 0.5, linewidth = 0.5,
                          color = "grey")
        ax = sns.lineplot(data = averageDF, x = "Generation", y = "MI_70 org frequency",
                          color = "blue")
        ax.set_ylabel("SI strategists' frequency")
    ax.grid(True)    
    return ax.get_figure()

def vectorField(resultsDir="./Output/gridResults/",csvfileName = "MIorg", equilibriumDir="./Output/gridEquilibrium/0.3MI/", startDir = "./Output/gridEquilibrium/0.7MI/",num_replicas = 10, axis=None, verbose=False):
    """
    Plotting function for vector field

    Parameters
    ----------
    resultsDir : str
        Path to grid results. Must follow same format as gridResults in repo.
    csvfileName : str
        File name without csv
    equilibriumDir : str
        Path to equilibrium runs. Should follow same format as repo
    num_replicas : int
        Num replicas for each grid

    Axis: axes matploblib object or None
        Used for main figures

    verbose : bool
        If True, also prints the start and end coordinates of every grid cell
    Returns
    -------
    figure
    """
    subFolders_AUPRC = ["0.1AUPRC","0.3AUPRC","0.5AUPRC","0.7AUPRC","0.9AUPRC"]
    
    if axis is None:
        fig, ax = plt.subplots()
    else:
        fig = plt.gcf()
        ax = axis
    
    #Computes equilbrium point
    equilbriumDFs = [pd.read_csv(f"{equilibriumDir}/{csvfileName}_{i}.csv") for i in range(num_replicas)]
    terminalYvals = [df["Mean AUPRC"].iloc[-1] for df in equilbriumDFs]
    terminalXvals = [100*(1 - df["Mean BS %MI"].iloc[-1]) for df in equilbriumDFs]
    terminalCoordinates = (mean(terminalXvals), mean(terminalYvals))
    
    #Computes steady state w/o selection
    startDFs = [pd.read_csv(f"{startDir}/{csvfileName}_{i}.csv") for i in range(num_replicas)]
    startYvals = [df["Mean AUPRC"].iloc[-1] for df in startDFs]
    startXvals = [100*(1 - df["Mean BS %MI"].iloc[-1]) for df in startDFs]
    startCoordinates = (mean(startXvals), mean(startYvals))  
    print("start coordinates: ", startCoordinates)
    print("end coordinates: ", terminalCoordinates)
    
    for subfolder_AUPRC in subFolders_AUPRC:
        subFolders = ["0.1MI","0.3MI","0.5MI","0.7MI","0.9MI"]
        
        dfs = []
        for subfolder in subFolders:
            replicaDF = []
            for i in range(num_replicas):
                filePath = f"{resultsDir}/{subfolder_AUPRC}/{subfolder}/{csvfileName}_{i}.csv"
                if os.path.exists(filePath):
                    replicaDF.append(pd.read_csv(filePath))
            dfs.append(replicaDF)

        vectors = {}
        
        #
        for i in range(len(subFolders)):
            #Grabs starting and ending coordinates for each replica
            starting_position_avg = [0,0] #Average starting position (x,y)
            ending_position_avg = [0,0] #Average ending position(x,y)
            for k in range(len(dfs[i])):
                df = dfs[i][k]
                starting_position_x = df["Mean BS %MI"].iloc[0]
                starting_position_y = df["Mean AUPRC"].iloc[0]
                ending_position_x = df["Mean BS %MI"].iloc[-1]
                ending_position_y = df["Mean AUPRC"].iloc[-1]
                
                starting_position_x = 100*(1 - starting_position_x) #Changes from %MI to %PI
                ending_position_x = 100*(1 - ending_position_x)
                
                starting_position_avg[0] += starting_position_x
                starting_position_avg[1] += starting_position_y
                ending_position_avg[0] += ending_position_x
                ending_position_avg[1] += ending_position_y
            
            #Normalizes 
            starting_position_avg = [coord/num_replicas for coord in starting_position_avg]
            ending_position_avg = [coord/num_replicas for coord in ending_position_avg]
            vectors[subFolders[i]] = (starting_position_avg, ending_position_avg)
    
            
        if verbose:
            print(subfolder_AUPRC)

        for key in vectors.keys():
            ax.plot([vectors[key][0][0]], [vectors[key][0][1]], 'ro-')
            #ax.plot([vectors[key][i][0] for i in range(2)], 'ro-')
            if verbose:
                print(key)
                for i in range(2):
                    print(f"({vectors[key][i][0]},{vectors[key][i][1]})")
            ax.annotate(
                '', 
                xy=(vectors[key][1][0], vectors[key][1][1]), xytext=(vectors[key][0][0], vectors[key][0][1]),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5,shrinkA=0, shrinkB=0), zorder=2
            )
    ax.scatter(startCoordinates[0], startCoordinates[1], marker='x', s=150, c='black', linewidths=2, label='Special Point', zorder=5,clip_on=False)
    ax.scatter(terminalCoordinates[0], terminalCoordinates[1], marker='*', s=300, c='black', edgecolors='white', linewidths=0.3, label='Special Point',zorder=5)
    ax.set_xlim(0, 100) 
    ax.set_ylim(0, 1) 
    ax.set_xlabel("Percent SI")
    ax.set_ylabel("AUPRC")
    ax.grid(True)

    return fig

def evolvabilityFig(
    subdirs=("./Output/ICreevolution/", "./Output/MIreevolution/"),
    target_len=1000,
    ax=None,
    ymin=0.0,
    shade=True
):
    """
    Plot Top AUPRC vs Generation for IC and MI re-evolution runs.

    Generates:
        (1) A combined plot (IC=red, MI=blue)
        (2) A second figure with vertical subplots (IC on top, MI on bottom).

    Parameters
    ----------
    subdirs : (str, str)
        Tuple of (IC_root, MI_root) directories to search recursively for CSVs.
    target_len : int
        Number of generations to display. Each run is padded with its last value
        (or truncated) to this length.
    ax : matplotlib.axes.Axes or None
        If provided, draw combined plot on this axes; otherwise create one.
    ymin : float
        Minimum y-axis value (default 0.0).
    shade : bool
        Whether to include IQR shading (default True).

    Returns
    -------
    tuple:
        (ax_combined, fig_separate)
    """

    IC_root, MI_root = subdirs

    # --- find CSVs ---
    def find_csvs(root):
        paths = []
        for r, _, files in os.walk(root):
            for f in files:
                if f.lower().endswith(".csv"):
                    paths.append(os.path.join(r, f))
        return sorted(paths)

    # --- extract "Top AUPRC" and pad/truncate ---
    def series_from_csv(path):
        df = pd.read_csv(path)
        if "Top AUPRC" not in df.columns:
            raise KeyError(f"'Top AUPRC' column not found in {path}")
        y = df["Top AUPRC"].to_numpy(dtype=float)
        if len(y) == 0:
            y = np.zeros(target_len)
        elif len(y) < target_len:
            y = np.concatenate([y, np.full(target_len - len(y), y[-1])])
        else:
            y = y[:target_len]
        return y

    # --- load data ---
    IC_runs = [series_from_csv(p) for p in find_csvs(IC_root)]
    MI_runs = [series_from_csv(p) for p in find_csvs(MI_root)]
    if len(IC_runs) == 0 and len(MI_runs) == 0:
        raise RuntimeError("No CSV files found under the provided subdirectories.")

    x = np.arange(target_len)

    # --- helper for plotting one group ---
    def plot_bundle(ax, runs, color, label_median):
        if len(runs) == 0:
            return
        Y = np.vstack(runs)

        # faint individual runs
        for y in Y:
            ax.plot(x, y, color=color, alpha=0.15, linewidth=0.8)

        # median & IQR
        median = np.nanmedian(Y, axis=0)
        q25 = np.nanpercentile(Y, 25, axis=0)
        q75 = np.nanpercentile(Y, 75, axis=0)

        ax.plot(x, median, color=color, linewidth=2.5, label=label_median)
        if shade:
            ax.fill_between(x, q25, q75, color=color, alpha=0.2, linewidth=0)

    # --- Combined plot ---
    created_fig = False
    if ax is None:
        fig_combined, ax_combined = plt.subplots(figsize=(7, 4.5), dpi=120)
        created_fig = True
    else:
        ax_combined = ax

    plot_bundle(ax_combined, IC_runs, color="red",  label_median="IC median ± IQR")
    plot_bundle(ax_combined, MI_runs, color="blue", label_median="MI median ± IQR")

    ax_combined.set_xlabel("Generation")
    ax_combined.set_ylabel("Top AUPRC")
    ax_combined.set_xlim(0, target_len - 1)
    ax_combined.set_ylim(ymin, 1.0)
    ax_combined.legend(frameon=False)
    ax_combined.grid(True, alpha=0.2, linestyle="--")

    if created_fig:
        plt.tight_layout()

    # --- Separate (vertical) subplots ---
    fig_sep, (ax_ic, ax_mi) = plt.subplots(2, 1, figsize=(7, 8), dpi=120, sharex=True, sharey=True)

    # IC subplot (top)
    plot_bundle(ax_ic, IC_runs, color="red",  label_median="IC median ± IQR")
    ax_ic.set_title("IC Re-evolution")
    ax_ic.set_ylabel("Top AUPRC")
    ax_ic.set_xlim(0, target_len - 1)
    ax_ic.set_ylim(ymin, 1.0)
    ax_ic.grid(True, alpha=0.2, linestyle="--")

    # MI subplot (bottom)
    plot_bundle(ax_mi, MI_runs, color="blue", label_median="MI median ± IQR")
    ax_mi.set_title("MI Re-evolution")
    ax_mi.set_xlabel("Generation")
    ax_mi.set_ylabel("Top AUPRC")
    ax_mi.set_xlim(0, target_len - 1)
    ax_mi.set_ylim(ymin, 1.0)
    ax_mi.grid(True, alpha=0.2, linestyle="--")

    fig_sep.tight_layout()

    return ax_combined, fig_sep

def _get_ax(plt_ax):
    # If already an Axes, use it
    if hasattr(plt_ax, "plot"):
        return plt_ax

    # If it's an int, map to current figure axes
    fig = plt.gcf()
    axs = fig.axes

    if isinstance(plt_ax, int):
        return axs[plt_ax]

    # 
    return plt.gca()

def gpmMetricBoxplot(results_df, metric, plt_ax=0, ylabel=None):
    '''
    Function that saves boxplots (with scatter points) comparing motifs of
    ascending IC in terms of the specified `metric`.
    '''
    ax = _get_ax(plt_ax)
    
    # Group results by 'class' (percent IC)
    grouped = results_df.groupby('class')
    names, vals, xs = [], [] ,[]
    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf[metric].tolist())
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))
    
    # Enforce that results are sorted by ascending IC
    ic_percentages = [float(name) for name in names]
    sorted_indexes = np.argsort(ic_percentages)
    names_sorted = [names[i] for i in sorted_indexes]
    vals_sorted = [vals[i] for i in sorted_indexes]
    
    # Boxplots
    ax.boxplot(vals_sorted, tick_labels=names_sorted)
    ngroup = len(vals_sorted)
    clevels = np.linspace(0., 1., ngroup)
    
    # Scatter points
    for x, val, clevel in zip(xs, vals_sorted, clevels):
        ax.scatter(x, val, color=cm.prism(clevel), alpha=0.4)
    
    # Label for Y-axis
    if not ylabel:
        ylabel = metric
    ax.set_ylabel(ylabel)


def gpmOptimaHistogram(results_df, plt_ax=0, alignment='center', vlines=True):
    '''
    Function that saves vertical histograms comparing motifs of
    ascending IC in terms of the "NO" metric (Number of Optima).
    '''
    ax = _get_ax(plt_ax)
    if alignment not in ['center', 'left']:
        raise ValueError("alignment must be 'center' or 'left'.")
    
    grouped = results_df.groupby('class')
    names, vals = [], []
    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf['NO'].tolist())
    
    # Enforce that results are sorted by ascending IC
    ic_percentages = [float(name) for name in names]
    sorted_indexes = np.argsort(ic_percentages)
    names_sorted = [names[i] for i in sorted_indexes]
    vals_sorted = [vals[i] for i in sorted_indexes]
    
    # Computed quantities to aid plotting
    min_no, max_no = np.min(vals_sorted), np.max(vals_sorted)
    number_of_bins = max_no - min_no + 1
    binned_data_sets = [
        np.histogram(d, range=(min_no, max_no), bins=number_of_bins)[0]
        for d in vals_sorted
    ]
    # Define x locations
    # Older code
    '''
    binned_maximums = np.max(binned_data_sets, axis=1)
    # Minimum distance between categories to avoid overlap
    width = np.max(binned_maximums)
    # Add extra space to guarantee there is always some space in between
    extra_sp = width * 0.1
    width += extra_sp
    x_locations = np.arange(0, sum(binned_maximums)+width, width)
    '''
    # New code
    extra_sp = 10
    width = 100 + extra_sp
    x_locations = np.arange(0, 5*width, width)
    
    # The bin_edges are the same for all of the histograms
    bin_edges = np.linspace(min_no-0.5, max_no+0.5, number_of_bins + 1)
    heights = np.diff(bin_edges)
    #heights = heights * 0.98  # Add extra space between bars
    centers = bin_edges[:-1] + heights / 2
    
    ngroup = len(x_locations)
    clevels = np.linspace(0., 1., ngroup)
    
    # Cycle through and plot each histogram
    for x_loc, binned_data, clevel in zip(x_locations, binned_data_sets, clevels):
        if alignment == 'center':
            lefts = x_loc - 0.5 * binned_data
        elif alignment == 'left':
            lefts = x_loc
            if vlines:
                ax.vlines(x_loc, bin_edges[0], bin_edges[-1], linewidth=0.5, colors ='k')
        #lefts = x_loc - 0.5 * width
        #lefts = x_loc
        ax.barh(centers, binned_data, height=heights, left=lefts, color=cm.prism(clevel))
    
    # Plot averages
    averages = [np.mean(vals) for vals in vals_sorted]
    ax.plot(x_locations, averages, linestyle='dashed', color='k', linewidth=0.8)
    ax.scatter(x_locations, averages, color='k', marker='d')
    
    # Ticks and labels
    ax.set_xticks(x_locations, names_sorted)
    #axs[plt_ax].set_yticks(list(range(min_no, max_no+1)))
    ax.set_ylabel('number of optima')
    
    # Set x-axis limits
    if alignment == 'center':
        ax.set_xlim((x_locations[0]-0.5*width, x_locations[-1]+0.5*width))
    elif alignment == 'left':
       ax.set_xlim((x_locations[0]-extra_sp, x_locations[-1]+1.0*width))

def percentSIheatmapFig(axis=None, L=6, sample_size=1000, only_adj_MI=True,
                                   star_coords=(6.47, 75.92),
                                   cross_coords=(1.64, 26.51)):
    # 0. Create an axis if one was not supplied (standalone use)
    if axis is None:
        fig, axis = plt.subplots(figsize=(5.5, 4))

    # 1. Data Generation
    mvp = MI_VS_SI(L)
    num_Lmers = 4**L
    n_vals = [int(num_Lmers/(2**i)) for i in range(1, 2*L + 1)]
    
    df_list = []
    for n in n_vals:
        df = mvp.percent_SI_distribution_sample(n, sample_size=sample_size,
                                                only_adj_MI=only_adj_MI)
        df_list.append(df)
    final_df = pd.concat(df_list)

    # 2. Setup Bins
    TI_vals = sorted(list(set(int(x) for x in final_df['TI'])))
    x_bins = np.arange(0.5, max(TI_vals) + 1.5, 1)
    y_bins = np.linspace(0, 100, 9)

    # 3. Calculate Normalized Histogram
    h, xedges, yedges = np.histogram2d(final_df['TI'], final_df['percent_SI'],
                                       bins=[x_bins, y_bins])
    h_norm = h / (h.sum(axis=1)[:, None] + 1e-10)

    # 4. Plotting on the provided axis
    X, Y = np.meshgrid(xedges, yedges)
    im = axis.pcolormesh(X, Y, h_norm.T, cmap='viridis', vmin=0, vmax=1)

    # 5. Overlay Markers
    axis.scatter(star_coords[0], star_coords[1], marker='*', s=300, 
               c='black', edgecolors='white', linewidths=0.3, zorder=3)
    axis.scatter(cross_coords[0], cross_coords[1], marker='x', s=150, 
               c='black', edgecolors='white' ,linewidths=2, zorder=3)

    # 6. Cell Highlighting
    for coords in [star_coords, cross_coords]:
        i = np.searchsorted(xedges, coords[0], side='right') - 1
        j = np.searchsorted(yedges, coords[1], side='right') - 1
        if 0 <= i < len(xedges)-1 and 0 <= j < len(yedges)-1:
            axis.add_patch(patches.Rectangle(
                (xedges[i], yedges[0]), xedges[i+1]-xedges[i], yedges[-1]-yedges[0],
                linewidth=0.3, edgecolor='black', facecolor='none'))
            axis.add_patch(patches.Rectangle(
                (xedges[i], yedges[j]), xedges[i+1]-xedges[i], yedges[j+1]-yedges[j],
                linewidth=1.0, edgecolor='white', facecolor='none'))

    # 7. Labels and Aesthetics
    axis.set_xlabel('TI (bits)')
    axis.set_ylabel('Percent SI')
    axis.set_xticks(TI_vals)
    axis.set_title(f'L = {L}')
    
    # Retrieve the figure to add the colorbar
    fig = axis.get_figure()
    cbar = fig.colorbar(im, ax=axis)
    cbar.set_label('p(Percent SI | TI)', rotation=270, labelpad=20)

    return fig

def percentSIsweepFig(only_adj_MI=False, sample_size=1000, annotateGeneticCode=True):
    """
    Plotting function returning the multipane heatmap of the percent PI distribution
    across motif lengths (L from 2 to 8 bp)

    only_adj_MI: bool
        If True, percent PI is computed from PI and the MI between adjacent positions
        only. If False, the total correlation C is used in place of the adjacent MI

    sample_size: int
        Number of motifs sampled per TI value to estimate each distribution

    annotateGeneticCode: bool
        If True, marks the standard genetic code with a red diamond on the L = 3 panel

    Returns
    -------
    figureplot
    """
    make_fig_4C_S7_S8 = importlib.import_module("src.make_fig_4C_S7-8")

    L_values = list(range(2, 9))
    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    axes = axes.flatten()
    img = None

    #Last axis is reserved for the shared colorbar, so it gets no heatmap
    for ax, L in zip(axes[:-1], L_values):
        print(f"Generating data for L={L} ...")
        df = make_fig_4C_S7_S8.generate_data(L=L, only_adj_MI=only_adj_MI,
                                             use_sampling=True, sample_size=sample_size)

        annotations = None
        if L == 3 and annotateGeneticCode:
            annotations = make_fig_4C_S7_S8.FIG_S8_L3_ANNOTATIONS

        img = make_fig_4C_S7_S8.plot_percent_SI_heatmap(df, L=L, ax=ax, annotations=annotations,
                                                        vmin=0, vmax=sample_size)

    make_fig_4C_S7_S8.add_colorbar_in_empty_ax(fig, axes[-1], img, sample_size)
    fig.tight_layout()

    return fig

def gpmAnalysisFig(fileList=["./Output/GPM/2bp_motifs/GPM_analysis_2bp.csv",
                             "./Output/GPM/6bp_motifs/GPM_analysis_6bp.csv"]):
    """
    Plotting function returning the genotype-phenotype map analysis of motifs of
    ascending percent PI, measured through FDC, GDPDC and number of optima

    fileList: list
        2 element list with paths to GPM_analysis csv files
        Note: 2 column path must be specified first and 6 column path must be specified second

    Returns
    -------
    figureplot
    """
    df_2bp = pd.read_csv(fileList[0])
    df_6bp = pd.read_csv(fileList[1])

    fig, axs = plt.subplots(3, 2, figsize=(12, 7), sharey='row')

    gpmMetricBoxplot(df_2bp, 'FDC_r', plt_ax=axs[0, 0], ylabel='FDC')
    gpmMetricBoxplot(df_2bp, 'GDPDC_r', plt_ax=axs[1, 0], ylabel='GDPDC')
    gpmOptimaHistogram(df_2bp, plt_ax=axs[2, 0])

    gpmMetricBoxplot(df_6bp, 'FDC_r', plt_ax=axs[0, 1], ylabel='FDC')
    gpmMetricBoxplot(df_6bp, 'GDPDC_r', plt_ax=axs[1, 1], ylabel='GDPDC')
    gpmOptimaHistogram(df_6bp, plt_ax=axs[2, 1])

    axs[0, 0].set_title("2 column motifs")
    axs[0, 1].set_title("6 column motifs")

    #Shared y axis drops the right hand labels, but the ticks are kept for readability
    for ax in axs[:, 1]:
        ax.set_ylabel("")
        ax.tick_params(labelleft=True)

    axs[2, 0].set_xlabel("percent PI")
    axs[2, 1].set_xlabel("percent PI")
    fig.align_ylabels()

    return fig

def geneticCodeRobustnessFig(n_random=1000):
    """
    Plotting function returning the scatterplot of percent PI against robustness for
    the standard genetic code and a sample of randomized codes

    n_random: int
        Number of randomized genetic codes to sample

    Returns
    -------
    figureplot
    """
    make_fig_6B_C = importlib.import_module("src.make_fig_6B-C")

    return make_fig_6B_C.make_figure_6B(n_random=n_random, outputDir=None)

def geneticCodeLogosFig():
    """
    Plotting function returning the sequence logos of the 21 messages of the genetic
    code, with a bar reporting the total correlation C of each message

    Returns
    -------
    2 element tuple of figureplots, for the standard code and for a randomized code
    """
    make_fig_6B_C = importlib.import_module("src.make_fig_6B-C")

    return make_fig_6B_C.make_figure_6C()

def geneticCodeCubesFig():
    """
    Plotting function returning the sequence logos of the 21 messages of the standard
    genetic code, each paired with its genotype table drawn as a cube

    Returns
    -------
    figureplot
    """
    make_fig_S9 = importlib.import_module("src.make_fig_S9")

    return make_fig_S9.make_figure_S9(outputDir=None)