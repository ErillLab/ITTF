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

from Objects.genetic_code import GeneticCode
from helper_functions import compute_seq_logo_matrix


# -------------------------------------------------------------------------
# Fig. 6B - Scatterplot: %SI vs Robustness in standard/random genetic codes
# -------------------------------------------------------------------------

GenCode = GeneticCode()

# Percent SI and Robustness of standard genetic code
GC_perc_SI = GenCode.si_percent()
GC_robustness = GenCode.robustness()

# Percent SI and Robustness of random codes
rnd_perc_SI_list = []
rnd_robustness_list = []
for i in range(1000):
    GenCode.randomize_code()
    rnd_perc_SI_list.append(GenCode.si_percent())
    rnd_robustness_list.append(GenCode.robustness())

# Generate scatterplot
plt.figure(figsize=(4,4), dpi=600)
plt.scatter(rnd_robustness_list, rnd_perc_SI_list, alpha=0.2, label='Random codes')
plt.scatter([GC_robustness], [GC_perc_SI], marker='d', color='red', alpha=0.8, label='Genetic code')
plt.legend()
plt.ylabel('Percent SI')
plt.xlabel('Robustness')
min_y = min(50, min(rnd_perc_SI_list))
plt.ylim((min_y,100))
plt.xticks([15,20,25,30])

plt.savefig('genetic_code_percent_SI_VS_robustness.png', bbox_inches='tight')
plt.show()

plt.close()


# ----------------------------------------------------------------
# Fig. 6C - Sequence logos (3x7) for standard/random genetic codes
# ----------------------------------------------------------------


def plot_code_logos(genetic_code):
    
    # Extract data for the 21 messages (20 amino acids + Stop signal)
    data = {}
    for message in genetic_code.messages:
        codons = list(genetic_code.table[genetic_code.table.translations == message]['codons'])
        si, c = genetic_code._compute_info(codons)
        data[message] = (codons, c)
    # Sort amino acids alphabetically and leave 'Stop' signal at the end
    aa_list = list(set(data.keys()).difference(set(['Stop'])))
    aa_list.sort()
    sorted_messages = aa_list + ['Stop']
    
    # Plot logos
    fig, axes = plt.subplots(3, 7, figsize=(14, 7), sharey=True)
    for i, ax in enumerate(axes.flatten()):
        message = sorted_messages[i]
        seqs, C = data[message]  # Sequences and total correlation C
        df = compute_seq_logo_matrix(seqs)
        # Plot sequence logo for the three positions X, Y, Z
        logo = logomaker.Logo(df, ax=ax)
        # Add 4th bar (for total correlation C)
        ax.bar(3, C, width=0.6, color='darkslategrey')
        ax.axhline(2, color='black', linewidth=0.5, alpha=0.6)
        # Formatting
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
    plt.show()

# Plot standard genetic code
GenCode = GeneticCode()
plot_code_logos(GenCode)

# Plot one random Code
GenCode.randomize_code()
plot_code_logos(GenCode)

