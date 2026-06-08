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

sns.set_style("whitegrid")

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
    src.figureFunctions.Figure4solo(fileList2, axis=axes[1])

    # Remove individual legend from Figure4solo panel
    legend_b = axes[1].get_legend()
    if legend_b:
        legend_b.remove()

    axes[0].text(-0.1, 1.05, "A", transform=axes[0].transAxes,
             fontsize=16, fontweight="bold", va="top")

    axes[1].text(-0.1, 1.05, "B", transform=axes[1].transAxes,
                 fontsize=16, fontweight="bold", va="top")

    # Shared legend above panel B, top right
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='red', lw=2, label='SI'),
                       Line2D([0], [0], color='blue', lw=2, label='MI')]
    fig.legend(handles=legend_elements, title="Strategy", ncol=1,
               frameon=True, fontsize=11,
               loc='upper right', bbox_to_anchor=(0.98, 1.05))

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
    axs[2, 0].set_xlabel("Percent SI")
    axs[2, 1].set_xlabel("Percent SI")
    return fig
def MainFigure4_new(panelA_path = "./Figure5A.png",
                    inputDir = "./Output/GPM",
                    outputDir = "./Output/Figures/MainFigure4_new.pdf"):
    """
    Combined figure with panel A (image) on top and panel B (GPM analysis) below.

    Parameters
    ----------
    panelA_path : str
        Path to the PNG image for panel A
    inputDir : str
        Path to GPM results directory
    outputDir : str
        Output path for the saved PDF
    """
    import matplotlib.image as mpimg
    from matplotlib.gridspec import GridSpec

    df_2bp = pd.read_csv(inputDir + "/2bp_motifs/GPM_analysis_2bp.csv")
    df_6bp = pd.read_csv(inputDir + "/6bp_motifs/GPM_analysis_6bp.csv")

    img = mpimg.imread(panelA_path)

    # Normalize to [0,1] if needed
    if img.max() > 1.0:
        img = img.astype(np.float32) / 255.0

    # Crop leading whitespace using mean brightness per column (robust to single dark pixels)
    col_mean = img[:, :, :3].mean(axis=(0, 2))  # shape (W,)
    non_white = np.where(col_mean < 0.98)[0]
    left_crop = int(non_white[0]) if len(non_white) > 0 else 0
    print(f"[MainFigure4_new] left_crop={left_crop} / {img.shape[1]} px")
    img = img[:, left_crop:, :]

    fig_width = 12
    a_ratio = 4
    total_height = 18

    fig = plt.figure(figsize=(fig_width, total_height))
    gs = GridSpec(4, 2, figure=fig, height_ratios=[a_ratio, 1, 1, 1], hspace=0.15, wspace=0.15)

    # Panel A — image spanning full width
    ax_img = fig.add_subplot(gs[0, :])
    ax_img.imshow(img, aspect='auto')
    ax_img.axis("off")
    ax_img.text(-0.02, 1.05, "A", transform=ax_img.transAxes,
                fontsize=16, fontweight="bold", va="top")

    # Panel B subplots
    axs = np.array([[fig.add_subplot(gs[row, col]) for col in range(2)] for row in range(1, 4)])

    src.figureFunctions.generate_boxplots_new_mapanalysis(df_2bp, 'FDC_r',   plt_ax=axs[0, 0], ylabel='FDC')
    src.figureFunctions.generate_boxplots_new_mapanalysis(df_2bp, 'GDPDC_r', plt_ax=axs[1, 0], ylabel='GDPDC')
    src.figureFunctions.generate_histograms_for_NO_new_mapanalysis(df_2bp,   plt_ax=axs[2, 0])

    src.figureFunctions.generate_boxplots_new_mapanalysis(df_6bp, 'FDC_r',   plt_ax=axs[0, 1], ylabel='FDC')
    src.figureFunctions.generate_boxplots_new_mapanalysis(df_6bp, 'GDPDC_r', plt_ax=axs[1, 1], ylabel='GDPDC')
    src.figureFunctions.generate_histograms_for_NO_new_mapanalysis(df_6bp,   plt_ax=axs[2, 1])

    # Remove grid lines from all B subplots
    for ax in axs.flat:
        ax.grid(False)

    # Column titles
    axs[0, 0].set_title("L = 2", fontsize=14)
    axs[0, 1].set_title("L = 6", fontsize=14)

    # Panel B label on first plot row only (covers all three rows)
    axs[0, 0].text(-0.18, 1.08, "B", transform=axs[0, 0].transAxes,
                   fontsize=16, fontweight="bold", va="top")

    # Remove y-axis label on right column (keep ticks)
    for ax in axs[:, 1]:
        ax.set_ylabel("")
        ax.tick_params(labelleft=True)

    # Set GDPDC row y-axis limits with buffer, but only label -0.2, 0, 0.2
    axs[1, 0].set_ylim(-0.35, 0.35)
    axs[1, 0].set_yticks([-0.2, 0, 0.2])
    axs[1, 1].set_ylim(-0.35, 0.35)
    axs[1, 1].set_yticks([-0.2, 0, 0.2])

    axs[2, 0].set_xlabel("Percent SI")
    axs[2, 1].set_xlabel("Percent SI")

    # Align panel A so its left edge matches the "B" label position (which sits in the margin)
    fig.canvas.draw()
    ax0_pos = axs[0, 0].get_position()
    right_edge = axs[0, 1].get_position().x1
    # "B" label is at -0.18 axes coords relative to axs[0,0]
    left_edge = max(0, ax0_pos.x0 - 0.18 * ax0_pos.width)
    img_pos = ax_img.get_position()
    ax_img.set_position([left_edge, img_pos.y0, right_edge - left_edge, img_pos.height])

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