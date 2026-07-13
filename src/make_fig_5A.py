#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 04:46:02 2026

@author: elia
"""


import numpy as np
# import pandas as pd
# import itertools
import logomaker
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from helper_functions.helperFunctions import compute_seq_logo_matrix, total_information


# XXX Import compute_seq_logo_matrix from helperFunctions.py

# XXX Maybe move `total_information` to helperFunctions.py



# def compute_seq_logo_matrix(seqs):
#     N = len(seqs)
#     L = len(seqs[0])
#     alphabet = ['A', 'C', 'G', 'T']
#     mat = []
    
#     for pos in range(L):
#         counts = {nt: 0 for nt in alphabet}
#         for s in seqs:
#             counts[s[pos]] += 1
        
#         freqs = {nt: counts[nt] / N for nt in alphabet}
        
#         # Shannon entropy
#         H = -sum(p * np.log2(p) for p in freqs.values() if p > 0)
#         # Information (max entropy for a nucleotide is 2 bits)
#         R = 2 - H

#         # Height of each letter
#         h = {nt: freqs[nt] * R for nt in alphabet}
#         mat.append(h)
    
#     df = pd.DataFrame(mat)
#     return df


# def entr(counts):
#     N = sum(counts)
#     H = 0
#     for c in counts:
#         if c != 0:
#             H -= (c/N) * np.log2(c/N)
#     return H

# def all_l_mers(site):
#     if len(site) == 1:
#         return ['A','C','G','T']
#     elif len(site) == 2:
#         return ['AA','AC','AG','AT',
#                 'CA','CC','CG','CT',
#                 'GA','GC','GG','GT',
#                 'TA','TC','TG','TT']
#     else:
#         raise ValueError('!!!')

# def entropy(sites):
#     alph = all_l_mers(sites[0])
#     counts = [sites.count(symbol) for symbol in alph]
#     return entr(counts)

"""
def total_information(sites):
    ''' Returns the Total Information (TI = SI + C) for the input set of sequences. '''
    # A priori Shannon entropy for the L-mers
    L = len(sites[0])
    n_states = 4**L
    H_before = np.log2(n_states)
    # Observed Shannon entropy for the L-mers
    all_Lmers = [''.join(seq) for seq in itertools.product('TCAG', repeat=L)]
    counts = [sites.count(Lmer) for Lmer in all_Lmers]
    N = sum(counts)
    freqs = [c/N for c in counts]
    H_after = -sum(p * np.log2(p) for p in freqs if p > 0)
    return H_before - H_after

"""


BASE_Z = 0.0

def plot_affinity_landscape(ax, z, title):
    
    x = np.array([0, 1, 0, 1])
    y = np.array([0, 0, 1, 1])
    z = np.array(z)
    
    verts = np.column_stack([x, y, z])
    
    # Base square
    bx = [0, 1, 1, 0, 0]
    by = [0, 0, 1, 1, 0]
    bz = [BASE_Z] * 5
    
    ax.plot(bx, by, bz, color='grey', linewidth=2.5)
    
    # Vertical grey lines
    for xi, yi, zi in verts:
        ax.plot([xi, xi], [yi, yi], [BASE_Z, zi],
                color='grey', linewidth=2.5, zorder=2)
    
    # Surface
    tri1 = [verts[0], verts[2], verts[1]]
    tri2 = [verts[2], verts[1], verts[3]]
    
    if title == "No epistasis":
        l_alpha, r_alpha = 0.45, 0.45
    elif title == "Magnitude epistasis":
        l_alpha, r_alpha = 0.32, 0.55
    elif title == "Sign epistasis":
        l_alpha, r_alpha = 0.55, 0.35
    elif title == "Reciprocal sign epistasis":
        l_alpha, r_alpha = 0.70, 0.32
    else:
        raise ValueError('Unknown title!!!')
    ax.add_collection3d(Poly3DCollection([tri1], alpha=l_alpha, zorder=3))
    ax.add_collection3d(Poly3DCollection([tri2], alpha=r_alpha, zorder=3))
    
    # Base edges on top
    eps = 1e-3
    
    if title == "Magnitude epistasis":
        ax.plot([0,1], [0,0], [BASE_Z+eps]*2, color='grey', linewidth=2.5, zorder=10)
    
    if title in ["Sign epistasis", "Reciprocal sign epistasis"]:
        ax.plot([0,1], [0,0], [BASE_Z+eps]*2, color='grey', linewidth=2.5, zorder=10)
        ax.plot([1,1], [0,1], [BASE_Z+eps]*2, color='grey', linewidth=2.5, zorder=10)
        z_AC = z[1]
        ax.plot([1,1], [0,0], [BASE_Z+eps, z_AC], color='grey', linewidth=2.5, zorder=10)
    
    def draw_arrow(i, j):
        
        p0 = verts[i]
        p1 = verts[j]
        
        # Shaft
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]],
                color='black', linewidth=1.2, zorder=5)
        
        # Arrowhead
        mid = 0.5 * (p0 + p1)  # midpoint and direction
        v = p1 - p0
        
        # axis scaling compensation
        zmin, zmax = ax.get_zlim()
        # zscale = (zmax - zmin) / 1.0
        box_aspect = ax.get_box_aspect()
        zscale = (zmax - zmin) / box_aspect[2]
        
        scale_vec = np.array([1, 1, zscale])
        
        # scaled direction
        v_scaled = np.array([v[0], v[1], v[2] / zscale])
        v_scaled /= np.linalg.norm(v_scaled)
        
        # Arrowhead size (height and base)
        h = 0.12
        b = h / 4
        
        # choose surface triangle
        if {i, j}.issubset({0, 1, 2}):
            tri_pts = [verts[0], verts[2], verts[1]]
        else:
            tri_pts = [verts[2], verts[1], verts[3]]
        
        # Compute normal in scaled space
        t0, t1, t2 = tri_pts
        t0s = np.array([t0[0], t0[1], t0[2] / zscale])
        t1s = np.array([t1[0], t1[1], t1[2] / zscale])
        t2s = np.array([t2[0], t2[1], t2[2] / zscale])
        
        normal = np.cross(t1s - t0s, t2s - t0s)
        normal /= np.linalg.norm(normal)
        
        # Perpendicular direction in plane
        perp_scaled = np.cross(normal, v_scaled)
        perp_scaled /= np.linalg.norm(perp_scaled)
        
        # Construct triangle in scaled space
        mid_s = mid / scale_vec
        
        tip_s = mid_s + 0.5 * h * v_scaled
        base_center_s = mid_s - 0.5 * h * v_scaled
        
        p_base1_s = base_center_s + b * perp_scaled
        p_base2_s = base_center_s - b * perp_scaled
        
        # Map back to original space
        tip = tip_s * scale_vec
        p_base1 = p_base1_s * scale_vec
        p_base2 = p_base2_s * scale_vec
        
        # Draw triangle
        tri = [tip, p_base1, p_base2]
        
        ax.add_collection3d(Poly3DCollection([tri], color='black', zorder=6))
    
    # Draw edges
    draw_arrow(0,1)
    draw_arrow(0,2)
    draw_arrow(1,3)
    draw_arrow(2,3)
    
    # Axes
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis='z', labelsize=14)
    
    ax.set_box_aspect([1,1,0.7])
    ax.view_init(elev=25, azim=-60)
    
    # ax.set_proj_type('ortho')
    
    # Draw points
    ax.scatter(x, y, z, color='black', s=40, depthshade=False, zorder=12)
    labels = ['CC', 'AC', 'CA', 'AA']
    for xi, yi, zi, label in zip(x, y, z, labels):
        ax.text(xi, yi, zi + 0.015, label, fontsize=18, ha='center', zorder=10)
    
    # Title
    # ax.set_title(title)



# Create the four motifs (for the four types of epistasis)

NEm  = 42*['AA'] + 18*['AC'] + 28*['CA'] + 12*['CC']  # No epistasis
MEm  = 39*['AA'] + 35*['AC'] + 25*['CA'] +  1*['CC']  # Magnitude epistasis
SEm  = 48*['AA'] +  1*['AC'] + 26*['CA'] + 25*['CC']  # Sign epistasis
RSEm = 50*['AA'] +  4*['AC'] +  1*['CA'] + 45*['CC']  # Reciprocal sign epistasis
motifs_list = [NEm, MEm, SEm, RSEm]
titles = ["No epistasis", "Magnitude epistasis", "Sign epistasis", "Reciprocal sign epistasis"]

# Probabilities
#     'CC'  'AC'  'CA'  'AA'
z1 = [0.12, 0.18, 0.28, 0.42]
z2 = [0.01, 0.35, 0.25, 0.39]
z3 = [0.25, 0.01, 0.26, 0.48]
z4 = [0.45, 0.04, 0.01, 0.50]

# Affinity (PSSM) scores
z1 = [np.log2(p) + 4 for p in z1]
z2 = [np.log2(p) + 4 for p in z2]
z3 = [np.log2(p) + 4 for p in z3]
z4 = [np.log2(p) + 4 for p in z4]
z_list = [z1, z2, z3, z4]



# Plot SEQUENCE LOGOS with corresponding AFFINITY LANDSCAPES for the four types of epistasis

fig = plt.figure(figsize=(16, 8), dpi=500)

gs = GridSpec(2, 4, figure=fig,
              height_ratios=[0.5, 2],  # logo vs landscape
              width_ratios=[1, 1, 1, 1],
              hspace=0.03, wspace=0.4)  # height space, width space

for i in range(4):
    
    # Top row: SEQUENCE LOGOS
    ax_logo = fig.add_subplot(gs[0, i])
    seqs = motifs_list[i]
    df = compute_seq_logo_matrix(seqs)
    ax_logo.set_ylim(0, 2)
    logo = logomaker.Logo(df, ax=ax_logo)
    # Add Correlation bar to sequence logo
    TI = total_information(seqs)
    PI = df.sum().sum()
    C = TI - PI
    ax_logo.bar(2, C, color='darkslategrey')
    # Format sequence logo
    ax_logo.axhline(1, color='black', linewidth=0.5, alpha=0.6)
    ax_logo.axhline(2, color='black', linewidth=0.5, alpha=0.6)
    ax_logo.set_xlim(-0.5, 2.5)
    ax_logo.set_xticks([0, 1, 2])
    ax_logo.set_xticklabels(['X', 'Y', 'MI'])
    ax_logo.set_ylim(0, 2.2)
    ax_logo.set_yticks([0, 1, 2])
    ax_logo.tick_params(axis='both', labelsize=20)
    ax_logo.set_ylabel('Bits', fontsize=20)
    ax_logo.set_title(titles[i] + '\n', fontsize=14)
    
    # Bottom row: AFFINITY LANDSCAPES
    ax_land = fig.add_subplot(gs[1, i], projection='3d')
    plot_affinity_landscape(ax_land, z_list[i], titles[i])

#plt.tight_layout()
plt.savefig('XXX_test_epistasis_diagrams_merged.png', bbox_inches='tight')
plt.show()












