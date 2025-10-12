# -*- coding: utf-8 -*-
"""
Functions used to generate figures of publication
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
import src.figureFunctions

def MainFigure1(fileList1 = ["./Notebook/Results/Figure_1/informationDecay_ICpop_fig1.csv", "./Notebook/Results/Figure_1/informationDecay_MIpop_fig1.csv"],
                fileList2 = ["./Notebook/Results/Figure_3/informationDecay_ICbit.csv", "./Notebook/Results/Figure_3/informationDecay_MIbit.csv"]):
    """
    Function to generate main figure 1

    Parameters
    ----------
    fileList1 : List 
        List of paths to output csv files for 2bp information decay experiment
    fileList2 : List
        List of paths to output csv files for 6bp length information decay experiment
    Returns
    -------
    Figure
    """
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    
    #p1
    src.figureFunctions.Figure1(fileList1, axis=axes[0])
    
    #p2
    src.figureFunctions.Figure3(fileList2,axis=axes[1])
    
    return fig

def MainFigure2(fileList1 = ["./Notebook/Results/SuppFigure_1/informationDecay_mixedpop_1_fig1.csv", "./Notebook/Results/SuppFigure_1/informationDecay_mixedpop_0_fig1.csv"],
                fileList2 = ["./Notebook/Results/Figure_4/genomeSearch_6bitIC_12col.csv", "./Notebook/Results/Figure_4/genomeSearch_6bitMI_12col.csv"]):
    """
    Function to generate main figure 1

    Parameters
    ----------
    fileList1 : List 
        List of paths to output csv files for 2bp information decay experiment
    fileList2 : List
        List of paths to output csv files for 12bp genome search experiment 
    Returns
    -------
    Figure
    """
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    
    #p1
    src.figureFunctions.SupplementaryFigure2(fileList1, axis=axes[0])
    
    #p2
    src.figureFunctions.Figure4solo(fileList2,axis=axes[1])
    
    return fig    

def MainFigure3():
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    
    #p1
    src.figureFunctions.competitionFig(axis=axes[0])
    
    #p2
    src.figureFunctions.vectorField(axis=axes[1])
    
    return fig    

def SupplementaryFig1(fileList = [f"./Notebook/Results/SuppFigure_1/informationDecay_mixedpop_{percentageIC}_fig1.json" for percentageIC in [0,0.25,0.50,0.75,1]]):
    return src.figureFunctions.SupplementaryFigure1(fileList)

def SupplementaryFig2(fileList = [[f"./Notebook/Results/Figure_2/informationDecay_{TI}bitIC.csv", f"./Notebook/Results/Figure_2/informationDecay_{TI}bitMI.csv"] for TI in [0.5,1,1.5,2]]):
    return src.figureFunctions.Figure2(fileList)

def SupplementaryFig3(fileList = [[f"./Notebook/Results/Figure_4/genomeSearch_{int(L/2)}bitIC_{L}col.csv", f"./Notebook/Results/Figure_4/genomeSearch_{int(L/2)}bitMI_{L}col.csv"] for L in [10,12,14,16]]):
    return src.figureFunctions.Figure4(fileList)

def SupplementaryFig4(fileList = [[f"./Notebook/Results/Figure_5/genomeSearch_{int(information)}bitIC_12col.csv", f"./Notebook/Results/Figure_5/genomeSearch_{int(information)}bitMI_12col.csv"] for information in [5,6,7]]):
    return src.figureFunctions.Figure5(fileList)