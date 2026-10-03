"""

Code to generate Fig. 6B, 6C, 7B, and S9.

- Fig. 6B:
  Standard/random genetic code scatterplot of %SI vs Robustness

- Fig. 6C:
  Sequence logos-based representation of the standard genetic code and one random code.

- Fig. 7B:
  Sequence logos-based representation of the three 'idealized' motifs used for the
  analytical results, i.e., motifs SI(A), SI(B) and MI.

- Fig. S9:
  XXX

"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import logomaker
import warnings

from Objects.genetic_code import GeneticCode
from helper_functions.helperFunctions import compute_seq_logo_matrix


# -------------------------------------------------------------------------
# Fig. 6B - Scatterplot: %SI vs Robustness in standard/random genetic codes
# -------------------------------------------------------------------------

GenCode = GeneticCode()

# Percent SI and Robustness of standard genetic code
GC_perc_SI, perc_SI_dict = GenCode.si_percent()
GC_robustness = GenCode.robustness()

# Percent SI and Robustness of random codes
rnd_perc_SI_list = []
rnd_robustness_list = []
for i in range(1000):
    GenCode.randomize_code()
    rnd_perc_SI, _ = GenCode.si_percent()
    rnd_robustness = GenCode.robustness()
    rnd_perc_SI_list.append(rnd_perc_SI)
    rnd_robustness_list.append(rnd_robustness)


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


# --------------------------------------------------------
# Fig. S9 - Genetic code sequence logos with genetic cubes
# --------------------------------------------------------

BASES = ['U', 'C', 'A', 'G']

def normalize_codon(codon):
    ''' Uppercase a codon and map T -> U, so DNA-style triplets (e.g. 'GCT')
    line up with the RNA-style axis labels used by the cube ('GCU'). '''
    return codon.strip().upper().replace('T', 'U')

def cube_faces(origin, size=1.0):
    '''' Return the 6 quad faces of a unit cube placed at `origin`. '''
    x, y, z = origin
    s = size
    p = [(x, y, z), (x + s, y, z), (x + s, y + s, z), (x, y + s, z),
         (x, y, z + s), (x + s, y, z + s), (x + s, y + s, z + s), (x, y + s, z + s)]
    faces = [[p[0], p[1], p[2], p[3]],  # bottom
             [p[4], p[5], p[6], p[7]],  # top
             [p[0], p[1], p[5], p[4]],  # front
             [p[2], p[3], p[7], p[6]],  # back
             [p[1], p[2], p[6], p[5]],  # right
             [p[0], p[3], p[7], p[4]]]   # left
    return faces

def plot_amino_acid_cube(ax, codons, highlight_color='#ffc000',
                         elev=22, azim=-60, background_alpha=0.02, highlight_alpha=0.95,
                         label_axes=False, tick_fontsize=6, label_fontsize=9, zoom=0.82):
    '''
    Draw the 4x4x4 genetic code cube on an existing 3D Axes. Every cell is a
    true 1x1x1 cube, fully transparent except the cells listed in `codons`,
    which are drawn opaque in `highlight_color`.
    
    `codons`: iterable of codon triplets.
    
    `zoom`: shrinks the rendered cube (and its tick labels) inward within
    its axes box (<1 = smaller, more margin).
    
    Axes:
        x -> 1st codon position
        y -> 2nd codon position
        z -> 3rd codon position
    '''
    target_codons = {normalize_codon(c) for c in codons}
    
    for i1, b1 in enumerate(BASES):
        for i2, b2 in enumerate(BASES):
            for i3, b3 in enumerate(BASES):
                codon = b1 + b2 + b3
                faces = cube_faces((i1, i2, i3), size=1.0)
                
                if codon in target_codons:
                    facecolor = highlight_color
                    alpha = highlight_alpha
                    edgecolor = 'black'
                    linewidth = 0.5
                else:
                    facecolor = 'lightgray'
                    alpha = background_alpha
                    edgecolor = (0, 0, 0, min(1.0, background_alpha * 3))
                    linewidth = 0.15
                
                poly = Poly3DCollection(faces, facecolor=facecolor, edgecolor=edgecolor,
                                        linewidths=linewidth, alpha=alpha)
                ax.add_collection3d(poly)
    
    # True cubes: orthographic projection removes perspective foreshortening,
    # and equal box aspect + equal data limits keep all three axes at the
    # same scale.
    ax.set_proj_type('ortho')
    ax.set_xlim(0, 4); ax.set_ylim(0, 4); ax.set_zlim(0, 4)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    ax.view_init(elev=elev, azim=azim)
    
    ticks = [i + 0.5 for i in range(4)]
    ax.set_xticks(ticks); ax.set_xticklabels(BASES, fontsize=tick_fontsize)
    ax.set_yticks(ticks); ax.set_yticklabels(BASES, fontsize=tick_fontsize)
    ax.set_zticks(ticks); ax.set_zticklabels(BASES, fontsize=tick_fontsize)
    ax.invert_zaxis()
    
    if label_axes:
        ax.set_xlabel('1st base', fontsize=label_fontsize, labelpad=6)
        ax.set_ylabel('2nd base', fontsize=label_fontsize, labelpad=6)
        ax.set_zlabel('3rd base', fontsize=label_fontsize, labelpad=2)
    
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1, 1, 1, 0))
        axis.pane.set_edgecolor((0, 0, 0, 0.15))
    ax.grid(False)


def plot_code_logos_with_cubes(genetic_code, elev=22, azim=-60, highlight_color='#ffc000',
                               background_alpha=0.02, highlight_alpha=0.95, figsize=(20, 13),
                               logo_tick_fontsize=9, logo_label_fontsize=10,
                               title_fontsize=12, cube_tick_fontsize=6,
                               height_ratios=(1.3, 1.0), inner_hspace=0.35, cube_zoom=1):
    '''
    Generates a 3x7 multi-panel plot (21 panels: 20 amino acids + Stop).
    Each panel stacks a sequence logo (top) on top of that message's genetic cube (bottom).
    '''
    
    # Extract data for the 21 messages (20 amino acids + Stop signal)
    data = {}
    for message in genetic_code.messages:
        codons = list(genetic_code.table[genetic_code.table.translations == message]['codons'])
        si, c = genetic_code._compute_info(codons)
        data[message] = (codons, c)
    
    # Sort amino acids alphabetically and leave 'Stop' signal at the end
    aa_list = sorted(set(data.keys()).difference({'Stop'}))
    sorted_messages = aa_list + ['Stop']
    
    # Make plot
    fig = plt.figure(figsize=figsize)
    outer_grid = gridspec.GridSpec(3, 7, figure=fig, wspace=0.35, hspace=0.35)
    
    for i, message in enumerate(sorted_messages):
        seqs, C = data[message]  # Sequences and total correlation C
        df = compute_seq_logo_matrix(seqs)
        
        inner_grid = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer_grid[i], height_ratios=height_ratios, hspace=inner_hspace)
        ax_logo = fig.add_subplot(inner_grid[0])
        ax_cube = fig.add_subplot(inner_grid[1], projection='3d')
        
        # Make sure the logo panel is always in front of the cube panel
        ax_logo.set_zorder(ax_cube.get_zorder() + 1)
        ax_logo.patch.set_facecolor('white')
        ax_logo.patch.set_alpha(1.0)
        
        # Sequence logo (top)
        logomaker.Logo(df, ax=ax_logo)
        ax_logo.bar(3, C, width=0.6, color='darkslategrey')
        ax_logo.axhline(1, color='black', linewidth=0.5, alpha=0.6)
        ax_logo.tick_params(axis='x', labelsize=logo_tick_fontsize)
        ax_logo.tick_params(axis='y', labelsize=logo_tick_fontsize)
        ax_logo.set_xlim(-0.5, 3.5)
        ax_logo.set_xticks([0, 1, 2, 3])
        ax_logo.set_xticklabels(['X', 'Y', 'Z', 'C'])
        ax_logo.set_ylim(0, 2)
        if C > 3:
            warnings.warn(f'In the "{message}" logo the correlation bar exceeds the height '
                          f'of the plot: C={C}')
        ax_logo.set_ylabel('Bits' if i % 7 == 0 else '', fontsize=logo_label_fontsize)
        ax_logo.set_title(message, fontsize=title_fontsize)
        
        # Genetic cube (bottom)
        plot_amino_acid_cube(ax_cube, seqs, highlight_color=highlight_color,
                             elev=elev, azim=azim, background_alpha=background_alpha,
                             highlight_alpha=highlight_alpha, label_axes=False,
                             tick_fontsize=cube_tick_fontsize, zoom=cube_zoom)
        
        # Box (border) around messages that are encoded with inter-positional correlations
        if C > 0:
            pos_logo = ax_logo.get_position()
            pos_cube = ax_cube.get_position()
            x0 = min(pos_logo.x0, pos_cube.x0)
            x1 = max(pos_logo.x1, pos_cube.x1)
            y0 = min(pos_logo.y0, pos_cube.y0)
            y1 = max(pos_logo.y1, pos_cube.y1)
            # Enlarge box to include title and labels
            x0 -= 0.015
            x1 += 0.01
            y0 -= 0.02
            y1 += 0.021
            if i % 7 == 0:
                x0 -= 0.01
            rect = patches.Rectangle((x0, y0), x1-x0, y1-y0, transform=fig.transFigure, fill=False,
                                     edgecolor='darkslategrey', linewidth=3.5, zorder=1000)
            fig.add_artist(rect)
    return fig


GenCode = GeneticCode()
fig_S9 = plot_code_logos_with_cubes(GenCode, height_ratios=(0.8, 1.7))
fig_S9.savefig('fig_S9_genetic_code_logos_and_cubes.png', dpi=600, bbox_inches='tight')





