"""

Code to perform the Adaptive Landscape Analysis (ALA) of motifs.

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

from Objects.motif_FLA import MotifFLA
from Objects.motiforganism import MotifModel


# Input data
motifs_base_dirpath = '../data/MI_to_IC_2Col_Motifs'
ic_percent_classes = [0, 25, 50, 75, 100]  

# Output directory path
#results_dirpath = '../Output/FLA/'
results_dirpath = '../Output/ALA/'
os.makedirs(results_dirpath, exist_ok=True)


global_results = []

for ic_percent in ic_percent_classes:
    
    # "IC" here (for "Information content") is a synonim of PI (for "Positional Information")
    dirpath = '{}/IC_{}/'.format(motifs_base_dirpath, ic_percent)
    
    fdc_r_list = []
    fdc_pval_list = []
    gdpdc_r_list = []
    gdpdc_pval_list = []
    no_list = []
    for filename in os.listdir(dirpath):
        # Read Motif Model
        MM = MotifModel(None, None, 0.5, None, None, read_in=dirpath+filename)
        # Analyze Fitness Landscape
        FLA = MotifFLA(MM)
        fdc_r, fdc_pval = FLA.fdc()
        gdpdc_r, gdpdc_pval = FLA.gdpdc()
        no, o_kmers = FLA.no()
        # Compile lists
        fdc_r_list.append(fdc_r)
        fdc_pval_list.append(fdc_pval)
        gdpdc_r_list.append(gdpdc_r)
        gdpdc_pval_list.append(gdpdc_pval)
        no_list.append(no)
    
    # Class-specific results as a table
    df = pd.DataFrame({'FDC_r': fdc_r_list,
                       'FDC_pval': fdc_pval_list,
                       'GDPDC_r': gdpdc_r_list,
                       'GDPDC_pval': gdpdc_pval_list,
                       'NO': no_list})
    
    # Compile global results
    #class_name = 'PI_' + str(ic_percent)
    class_name = str(ic_percent)
    df['class'] = [class_name] * len(df)
    global_results.append(df)

# Results as table
results_df = pd.concat(global_results)
# Save results as CSV file
results_df.to_csv(results_dirpath + 'ALA_table.csv', index=False)

# --------------------------------------
# Plot results
# --------------------------------------


def generate_boxplots_new(results_df, metric, plt_ax=0, ylabel=None):
    '''
    Function that saves boxplots (with scatter points) comparing motifs of
    ascending IC in terms of the specified `metric`.
    '''
    # Group results by 'class' (percent IC)
    grouped = results_df.groupby('class')
    names, vals, xs = [], [] ,[]
    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf[metric].tolist())
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))
    
    # Enforce that results are sorted by ascending IC
    ic_percentages = [float(name) for name in names]
    sorted_indexes = np.argsort(ic_percentages)
    names_sorted = [names[i] for i in sorted_indexes]
    vals_sorted = [vals[i] for i in sorted_indexes]
    
    # Boxplots
    
    axs[plt_ax].boxplot(vals_sorted, tick_labels=names_sorted)
    ngroup = len(vals_sorted)
    clevels = np.linspace(0., 1., ngroup)
    
    # Scatter points
    for x, val, clevel in zip(xs, vals_sorted, clevels):
        axs[plt_ax].scatter(x, val, color=cm.prism(clevel), alpha=0.4)
    
    # Label for Y-axis
    if not ylabel:
        ylabel = metric
    axs[plt_ax].set_ylabel(ylabel)


def generate_histograms_for_NO_new(results_df, plt_ax=0, alignment='center', vlines=True):
    '''
    Function that saves vertical histograms comparing motifs of
    ascending IC in terms of the "NO" metric (Number of Optima).
    '''
    
    if alignment not in ['center', 'left']:
        raise ValueError("alignment must be 'center' or 'left'.")
    
    grouped = results_df.groupby('class')
    names, vals = [], []
    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf['NO'].tolist())
    
    # Enforce that results are sorted by ascending IC
    ic_percentages = [float(name) for name in names]
    sorted_indexes = np.argsort(ic_percentages)
    names_sorted = [names[i] for i in sorted_indexes]
    vals_sorted = [vals[i] for i in sorted_indexes]
    
    # Computed quantities to aid plotting
    min_no, max_no = np.min(vals_sorted), np.max(vals_sorted)
    number_of_bins = max_no - min_no + 1
    binned_data_sets = [
        np.histogram(d, range=(min_no, max_no), bins=number_of_bins)[0]
        for d in vals_sorted
    ]
    binned_maximums = np.max(binned_data_sets, axis=1)
    # Minimum distance between categories to avoid overlap
    width = np.max(binned_maximums)
    # Add extra space to guarantee there is always some space in between
    extra_sp = width * 0.1
    width += extra_sp
    x_locations = np.arange(0, sum(binned_maximums)+width, width)
    
    # The bin_edges are the same for all of the histograms
    bin_edges = np.linspace(min_no-0.5, max_no+0.5, number_of_bins + 1)
    heights = np.diff(bin_edges)
    #heights = heights * 0.98  # Add extra space between bars
    centers = bin_edges[:-1] + heights / 2
    
    ngroup = len(x_locations)
    clevels = np.linspace(0., 1., ngroup)
    
    # Cycle through and plot each histogram
    for x_loc, binned_data, clevel in zip(x_locations, binned_data_sets, clevels):
        if alignment == 'center':
            lefts = x_loc - 0.5 * binned_data
        elif alignment == 'left':
            lefts = x_loc
            if vlines:
                axs[plt_ax].vlines(x_loc, bin_edges[0], bin_edges[-1], linewidth=0.5, colors ='k')
        #lefts = x_loc - 0.5 * width
        #lefts = x_loc
        axs[plt_ax].barh(centers, binned_data, height=heights, left=lefts, color=cm.prism(clevel))
    
    # Plot averages
    averages = [np.mean(vals) for vals in vals_sorted]
    axs[plt_ax].plot(x_locations, averages, linestyle='dashed', color='k', linewidth=0.8)
    axs[plt_ax].scatter(x_locations, averages, color='k', marker='d')
    
    # Ticks and labels
    axs[plt_ax].set_xticks(x_locations, names_sorted)
    axs[plt_ax].set_yticks(list(range(min_no, max_no+1)))
    axs[plt_ax].set_ylabel('number of optima')
    
    # Set x-axis limits
    if alignment == 'center':
        axs[plt_ax].set_xlim((x_locations[0]-0.5*width, x_locations[-1]+0.5*width))
    elif alignment == 'left':
        axs[plt_ax].set_xlim((x_locations[0]-extra_sp, x_locations[-1]+1.0*width))



# -----------------------------------------------------------------------------
# Make Adaptive Landscape Analysis Figure (three panes)
# -----------------------------------------------------------------------------

fig, axs = plt.subplots(3, figsize=(6, 7))
generate_boxplots_new(results_df, 'FDC_r', plt_ax=0, ylabel='FDC')
generate_boxplots_new(results_df, 'GDPDC_r', plt_ax=1, ylabel='GDPDC')
generate_histograms_for_NO_new(results_df, plt_ax=2)
fig.align_ylabels()
plt.xlabel('percent PI')
# plt.show()
plt.savefig(results_dirpath + 'ALA_Figure.png', dpi=300)
plt.close()



# -----------------------------------------------------------------------------
# Old code to generate individual plots separately
# -----------------------------------------------------------------------------

"""

def generate_boxplots_old(results_df, metric, ylabel=None):
    '''
    Function that saves boxplots (with scatter points) comparing motifs of
    ascending IC in terms of the specified `metric`.
    '''
    # Group results by 'class' (percent IC)
    grouped = results_df.groupby('class')
    names, vals, xs = [], [] ,[]
    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf[metric].tolist())
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))
    
    # Enforce that results are sorted by ascending IC
    ic_percentages = [float(name.split('_')[1]) for name in names]
    sorted_indexes = np.argsort(ic_percentages)
    names_sorted = [names[i] for i in sorted_indexes]
    vals_sorted = [vals[i] for i in sorted_indexes]
    
    # Boxplots
    plt.boxplot(vals_sorted, tick_labels=names_sorted)
    ngroup = len(vals_sorted)
    clevels = np.linspace(0., 1., ngroup)
    
    # Scatter points
    for x, val, clevel in zip(xs, vals_sorted, clevels):
        plt.scatter(x, val, color=cm.prism(clevel), alpha=0.4)
    
    # Label for Y-axis
    if not ylabel:
        ylabel = metric
    plt.ylabel(ylabel)
    
    # Save figure
    plt.savefig(results_dirpath + metric + '.png', dpi=300)
    plt.close()

def generate_histograms_for_NO_old(results_df, alignment='center', vlines=True):
    '''
    Function that saves vertical histograms comparing motifs of
    ascending IC in terms of the "NO" metric (Number of Optima).
    '''
    
    if alignment not in ['center', 'left']:
        raise ValueError("alignment must be 'center' or 'left'.")
    
    grouped = results_df.groupby('class')
    names, vals = [], []
    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf['NO'].tolist())
    
    # Enforce that results are sorted by ascending IC
    ic_percentages = [float(name.split('_')[1]) for name in names]
    sorted_indexes = np.argsort(ic_percentages)
    names_sorted = [names[i] for i in sorted_indexes]
    vals_sorted = [vals[i] for i in sorted_indexes]
    
    # Computed quantities to aid plotting
    min_no, max_no = np.min(vals_sorted), np.max(vals_sorted)
    number_of_bins = max_no - min_no + 1
    binned_data_sets = [
        np.histogram(d, range=(min_no, max_no), bins=number_of_bins)[0]
        for d in vals_sorted
    ]
    binned_maximums = np.max(binned_data_sets, axis=1)
    # Minimum distance between categories to avoid overlap
    width = np.max(binned_maximums)
    # Add extra space to guarantee there is always some space in between
    extra_sp = width * 0.1
    width += extra_sp
    x_locations = np.arange(0, sum(binned_maximums)+width, width)
    
    # The bin_edges are the same for all of the histograms
    bin_edges = np.linspace(min_no-0.5, max_no+0.5, number_of_bins + 1)
    heights = np.diff(bin_edges)
    #heights = heights * 0.98  # Add extra space between bars
    centers = bin_edges[:-1] + heights / 2
    
    ngroup = len(x_locations)
    clevels = np.linspace(0., 1., ngroup)
    
    # Cycle through and plot each histogram
    fig, ax = plt.subplots()
    for x_loc, binned_data, clevel in zip(x_locations, binned_data_sets, clevels):
        if alignment == 'center':
            lefts = x_loc - 0.5 * binned_data
        elif alignment == 'left':
            lefts = x_loc
            if vlines:
                ax.vlines(x_loc, bin_edges[0], bin_edges[-1], linewidth=0.5, colors ='k')
        #lefts = x_loc - 0.5 * width
        #lefts = x_loc
        ax.barh(centers, binned_data, height=heights, left=lefts, color=cm.prism(clevel))
    
    # Plot averages
    averages = [np.mean(vals) for vals in vals_sorted]
    ax.plot(x_locations, averages, linestyle='dashed', color='k', linewidth=0.8)
    ax.scatter(x_locations, averages, color='k', marker='d')
    
    # Ticks and labels
    ax.set_xticks(x_locations, names_sorted)
    ax.set_yticks(list(range(min_no, max_no+1)))
    ax.set_ylabel('number of optima')
    
    # Set x-axis limits
    if alignment == 'center':
        ax.set_xlim((x_locations[0]-0.5*width, x_locations[-1]+0.5*width))
    elif alignment == 'left':
        ax.set_xlim((x_locations[0]-extra_sp, x_locations[-1]+1.0*width))
    
    # Save figure
    if alignment == 'left' and vlines:
        vlines_tag = '_vlines'
    else:
        vlines_tag = ''
    plt.savefig(results_dirpath + 'NO_hist_' + alignment + vlines_tag + '.png', dpi=300)
    plt.close()


generate_boxplots_old(results_df, 'GDPDC_r', ylabel='genotype distace - phenotype distance correlation')
generate_boxplots_old(results_df, 'FDC_r', ylabel='fitness - distance from optimum correlation')
generate_histograms_for_NO_old(results_df, alignment='center')

"""
















