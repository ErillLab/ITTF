"""

Code to perform the Fitness Landscape Analysis (FLA) of motifs.

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
results_dirpath = '../Output/FLA/'
os.makedirs(results_dirpath, exist_ok=True)


global_results = []

for ic_percent in ic_percent_classes:
    
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
    class_name = 'IC_' + str(ic_percent)
    df['class'] = [class_name] * len(df)
    global_results.append(df)

# Results as table
results_df = pd.concat(global_results)
# Save results as CSV file
results_df.to_csv(results_dirpath + 'FLA_table.csv', index=False)

# --------------------------------------
# Plot results
# --------------------------------------

def generate_boxplots(results_df, metric, ylabel=None):
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




generate_boxplots(results_df, 'GDPDC_r', ylabel='genotype distace - phenotype distance correlation')
generate_boxplots(results_df, 'FDC_r', ylabel='fitness - distance from optimum correlation')
generate_boxplots(results_df, 'NO', ylabel='number of optima')






























