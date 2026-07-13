# -*- coding: utf-8 -*-
"""

Code to generate Fig. 7B.

Fig. 7B: Sequence logos-based representation of the three 'idealized' motifs
used for the analytical results, i.e., motifs SI(A), SI(B) and MI.

"""

import matplotlib.pyplot as plt
import logomaker

from helper_functions.helperFunctions import compute_seq_logo_matrix


SIAm = ['AA', 'AC', 'AG', 'AT']
SIBm = ['AA', 'AC', 'CA', 'CC']
MIm = ['AC', 'CG', 'GT', 'TA']
motif_names = ['SI(A)', 'SI(B)', 'MI']
seqs_l = [SIAm, SIBm, MIm]
C_l = [0, 0, 2]

# Plot logos
fig, axes = plt.subplots(1, 3, figsize=(14, 2.5), sharey=False)
for i, ax in enumerate(axes.flatten()):
    seqs = seqs_l[i]
    df = compute_seq_logo_matrix(seqs)
    # Plot sequence logo for the three positions X, Y, Z
    ax.set_ylim(0, 2)
    logo = logomaker.Logo(df, ax=ax)
    # Add 4th bar (for total correlation C)
    ax.bar(2, C_l[i], width=0.6, color='darkslategrey')
    ax.axhline(2, color='black', linewidth=0.5, alpha=0.6)
    ax.axhline(1, color='black', linewidth=0.5, alpha=0.6)
    # Formatting
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_xlim(-0.5, 2.5)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['X', 'Y', 'MI'])
    ax.set_ylim(0, 2.2)
    ax.set_yticks([0, 1, 2])
    ax.set_ylabel('Bits', fontsize=22)
    #ax.set_title(motif_names[i], fontsize=22)
plt.tight_layout()
plt.show()







