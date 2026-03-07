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

def MainFigure2(fileList1 = ["./Notebook/Results/Figure_1/informationDecay_ICpop_fig1.csv", "./Notebook/Results/Figure_1/informationDecay_MIpop_fig1.csv"],
                fileList2 = ["./Notebook/Results/Figure_4/genomeSearch_6bitIC_12col.csv", "./Notebook/Results/Figure_4/genomeSearch_6bitMI_12col.csv"],
                outputDir = "./Output/Figures/MainFigure2.pdf"):
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
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    
    #p1
    src.figureFunctions.Figure1(fileList1, axis=axes[0])
    
    #p2
    src.figureFunctions.Figure4solo(fileList2,axis=axes[1])
    
    if outputDir:
        fig.savefig(outputDir, format="pdf", bbox_inches="tight")
    
    return fig

def MainFigure3(outputDir = "./Output/Figures/MainFigure3.pdf"):
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    
    #p1
    src.figureFunctions.competitionFig(axis=axes[0])
    
    #p2
    src.figureFunctions.vectorField(axis=axes[1])
    
    if outputDir:
        fig.savefig(outputDir, format="pdf", bbox_inches="tight")
        
    return fig    

def SupplementaryFig1(
    fileList=[f"./Notebook/Results/SuppFigure_1/informationDecay_mixedpop_{percentageIC}_fig1.csv"
              for percentageIC in [0, 0.25, 0.50, 0.75, 1]],
    outputDir="./Output/Figures/SuppFigure1.pdf"
):
    fig = src.figureFunctions.SupplementaryFigure1(fileList)
    fig.savefig(outputDir, format="pdf", bbox_inches="tight")
    return fig


def SupplementaryFig2(
    fileList=[[f"./Notebook/Results/Figure_2/informationDecay_{TI}bitIC.csv",
               f"./Notebook/Results/Figure_2/informationDecay_{TI}bitMI.csv"]
              for TI in [0.5, 1, 1.5, 2]],
    outputDir="./Output/Figures/SuppFigure2.pdf"
):
    fig = src.figureFunctions.Figure2(fileList)
    fig.savefig(outputDir, format="pdf", bbox_inches="tight")
    return fig


def SupplementaryFig3(
    fileList=["./Notebook/Results/Figure_3/informationDecay_ICbit.csv",
              "./Notebook/Results/Figure_3/informationDecay_MIbit.csv"],
    outputDir="./Output/Figures/SuppFigure3.pdf"
):
    fig = src.figureFunctions.Figure3(fileList)
    fig.savefig(outputDir, format="pdf", bbox_inches="tight")
    return fig


def SupplementaryFig4(
    fileList=[[f"./Notebook/Results/Figure_4/genomeSearch_{int(L/2)}bitIC_{L}col.csv",
               f"./Notebook/Results/Figure_4/genomeSearch_{int(L/2)}bitMI_{L}col.csv"]
              for L in [10, 12, 14, 16]],
    outputDir="./Output/Figures/SuppFigure4.pdf"
):
    fig = src.figureFunctions.Figure4(fileList)
    fig.savefig(outputDir, format="pdf", bbox_inches="tight")
    return fig


def SupplementaryFig5(
    fileList=[[f"./Notebook/Results/Figure_5/genomeSearch_{int(info)}bitIC_12col.csv",
               f"./Notebook/Results/Figure_5/genomeSearch_{int(info)}bitMI_12col.csv"]
              for info in [5, 6, 7]],
    outputDir="./Output/Figures/SuppFigure5.pdf"
):
    fig = src.figureFunctions.Figure5(fileList)
    fig.savefig(outputDir, format="pdf", bbox_inches="tight")
    return fig