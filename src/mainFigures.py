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
    
    axes[0].text(-0.1, 1.05, "A", transform=axes[0].transAxes,
             fontsize=16, fontweight="bold", va="top")

    axes[1].text(-0.1, 1.05, "B", transform=axes[1].transAxes,
                 fontsize=16, fontweight="bold", va="top")

    if outputDir:
        fig.savefig(outputDir, format="pdf", bbox_inches="tight")
    
    return fig

def MainFigure3(outputDir = "./Output/Figures/MainFigure3.pdf"):
    # Updated to 3 columns and wider width
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    
    # Panel 1
    src.figureFunctions.competitionFig(axis=axes[0])
    
    # Panel 2
    src.figureFunctions.vectorField(axis=axes[1])
    
    # Panel 3: The Heatmap
    # Note: Assuming MI_VS_PI and generate_representation_heatmap are accessible
    src.figureFunctions.generate_representation_heatmap(axis=axes[2])
    
    # Optional: Add Panel Labels (A, B, C) if you want consistency with Fig 4
    for i, label in enumerate(['A', 'B', 'C']):
        axes[i].text(-0.15, 1.1, label, transform=axes[i].transAxes,
                     fontsize=16, fontweight="bold", va="top")
    
    if outputDir:
        fig.savefig(outputDir, format="pdf", bbox_inches="tight")
        
    return fig
def MainFigure4(inputDir = "./Output/GPM"):
    df_2bp = pd.read_csv(inputDir + "/2bp_motifs/GPM_analysis_2bp.csv")
    df_6bp = pd.read_csv(inputDir + "/6bp_motifs/GPM_analysis_6bp.csv")
    fig, axs = plt.subplots(3, 2, figsize=(12, 7))
    
    #Left column (2bp)
    src.figureFunctions.generate_boxplots_new_mapanalysis(df_2bp, 'FDC_r', plt_ax=axs[0, 0], ylabel='FDC')
    src.figureFunctions.generate_boxplots_new_mapanalysis(df_2bp, 'GDPDC_r', plt_ax=axs[1, 0], ylabel='GDPDC')
    src.figureFunctions.generate_histograms_for_NO_new_mapanalysis(df_2bp, plt_ax=axs[2, 0])
    
    #right column (6bp)
    src.figureFunctions.generate_boxplots_new_mapanalysis(df_6bp, 'FDC_r', plt_ax=axs[0, 1], ylabel='FDC')
    src.figureFunctions.generate_boxplots_new_mapanalysis(df_6bp, 'GDPDC_r', plt_ax=axs[1, 1], ylabel='GDPDC')
    src.figureFunctions.generate_histograms_for_NO_new_mapanalysis(df_6bp, plt_ax=axs[2, 1])
    
    
    # Column titles
    axs[0, 0].set_title("L = 2")
    axs[0, 1].set_title("L = 6")
    
    # Panel labels for left column
    axs[0, 0].text(-0.18, 1.08, "A", transform=axs[0, 0].transAxes,
                   fontsize=16, fontweight="bold", va="top")
    axs[1, 0].text(-0.18, 1.08, "B", transform=axs[1, 0].transAxes,
                   fontsize=16, fontweight="bold", va="top")
    axs[2, 0].text(-0.18, 1.08, "C", transform=axs[2, 0].transAxes,
                   fontsize=16, fontweight="bold", va="top")
    
    # Remove y-axis LABEL only (keep ticks)
    for ax in axs[:, 1]:
        ax.set_ylabel("")
        ax.tick_params(labelleft=True)
    
    # Bottom x labels
    axs[2, 0].set_xlabel("percent PI")
    axs[2, 1].set_xlabel("percent PI")
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