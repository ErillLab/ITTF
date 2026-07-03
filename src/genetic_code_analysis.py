"""

Code to generate Fig. 6B, 6C and 7B.

- Fig. 6B:
  Standard/random genetic code scatterplot of %SI vs Robustness

- Fig. 6C:
  Sequence logos-based representation of the standard genetic code and one random code.

- Fig. 7B:
  Sequence logos-based representation of the three 'idealized' motifs used for the
  analytical results, i.e., motifs SI(A), SI(B) and MI.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker

from Objects.genetic_code import GeneticCode


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

def compute_information_matrix(seqs):
    N = len(seqs)
    alphabet = ['A', 'C', 'G', 'T']
    mat = []
    for pos in range(len(seqs[0])):
        counts = {nt: 0 for nt in alphabet}
        for s in seqs:
            counts[s[pos]] += 1
        freqs = {nt: counts[nt] / N for nt in alphabet}
        # Shannon entropy
        H = -sum(p * np.log2(p) for p in freqs.values() if p > 0)
        SI = 2 - H
        # Height of each letter at this position of the sequence logo
        mat.append({nt: freqs[nt] * SI for nt in alphabet})
    return pd.DataFrame(mat)


def plot_code_logos(genetic_code):
    
    # Extract data
    messages_l = []
    codons_l = []
    C_l = []
    for message in genetic_code.messages:
        messages_l.append(message)
        codons = list(genetic_code.table[genetic_code.table.translations == message]['codons'])
        codons_l.append(codons)
        si, c = genetic_code._compute_info(codons)
        C_l.append(c)
    
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
        df = compute_information_matrix(seqs)
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

# ------------------------------------------------------------------------------
# Fig. 7B - Sequence logos (1x3) for the three idealized cases: SI(A), SI(B), MI
# ------------------------------------------------------------------------------


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
    df = compute_information_matrix(seqs)
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


















