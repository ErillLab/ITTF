# -*- coding: utf-8 -*-
"""

Code to generate Fig. 6B and 6C

- Fig. 6B:
  Standard/random genetic code scatterplot of %SI vs Robustness

- Fig. 6C:
  Sequence logos-based representation of the standard genetic code and one random code.

"""

import matplotlib.pyplot as plt
import logomaker
import warnings

from src.Objects.genetic_code import GeneticCode
from src.helper_functions.helperFunctions import compute_seq_logo_matrix


# -------------------------------------------------------------------------
# Fig. 6B - Scatterplot: %SI vs Robustness in standard/random genetic codes
# -------------------------------------------------------------------------

def make_figure_6B(n_random=1000, outputDir='genetic_code_percent_SI_VS_robustness.png'):
    GenCode = GeneticCode()
    
    GC_perc_SI, perc_SI_dict = GenCode.si_percent()
    GC_robustness, robustness_dict = GenCode.robustness()
    
    rnd_perc_SI_list = []
    rnd_robustness_list = []
    for i in range(n_random):
        if i % 100 == 0:
            print(f"Figure 6B: {i}/{n_random} random codes...")
        GenCode.randomize_code()
        rnd_perc_SI, _ = GenCode.si_percent()
        rnd_robustness, _ = GenCode.robustness()
        rnd_perc_SI_list.append(rnd_perc_SI)
        rnd_robustness_list.append(rnd_robustness)
    
    fig, ax = plt.subplots(figsize=(4, 4), dpi=600)
    ax.scatter(rnd_robustness_list, rnd_perc_SI_list, alpha=0.2, label='Random codes')
    ax.scatter([GC_robustness], [GC_perc_SI], marker='d', color='red', alpha=0.8, label='Genetic code')
    ax.legend()
    ax.set_ylabel('Percent SI')
    ax.set_xlabel('Robustness')
    min_y = min(50, min(rnd_perc_SI_list))
    ax.set_ylim((min_y, 100))
    ax.set_xticks([15, 20, 25, 30])
    
    if outputDir:
        fig.savefig(outputDir, bbox_inches='tight')
    
    return fig


# ----------------------------------------------------------------
# Fig. 6C - Sequence logos (3x7) for standard/random genetic codes
# ----------------------------------------------------------------

def plot_code_logos(genetic_code):
    data = {}
    for message in genetic_code.messages:
        codons = list(genetic_code.table[genetic_code.table.translations == message]['codons'])
        si, c = genetic_code._compute_info(codons)
        data[message] = (codons, c)
    aa_list = list(set(data.keys()).difference(set(['Stop'])))
    aa_list.sort()
    sorted_messages = aa_list + ['Stop']
    
    fig, axes = plt.subplots(3, 7, figsize=(14, 7), sharey=True)
    for i, ax in enumerate(axes.flatten()):
        message = sorted_messages[i]
        seqs, C = data[message]
        df = compute_seq_logo_matrix(seqs)
        logo = logomaker.Logo(df, ax=ax)
        ax.bar(3, C, width=0.6, color='darkslategrey')
        ax.axhline(2, color='black', linewidth=0.5, alpha=0.6)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        ax.set_xlim(-0.5, 3.5)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['X', 'Y', 'Z', 'C'])
        ax.set_ylim(0, 3)
        if C > 3:
            warnings.warn(f'In the "{message}" logo the correlation bar exceeds the height '
                          f'of the plot: C={C}')
        ax.set_ylabel('Bits' if i % 7 == 0 else '', fontsize=22)
        ax.set_title(message, fontsize=22)
    plt.tight_layout()
    return fig


def make_figure_6C():
    GenCode = GeneticCode()
    fig_standard = plot_code_logos(GenCode)
    
    GenCode.randomize_code()
    fig_random = plot_code_logos(GenCode)
    
    return fig_standard, fig_random


def main():
    make_figure_6B()
    make_figure_6C()



if __name__ == "__main__":
    main()









