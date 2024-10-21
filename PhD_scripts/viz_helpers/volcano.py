import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
import matplotlib.patheffects as PathEffects
import os


def volcano(data, log2fc='log2FoldChange', pvalue='padj', symbol='symbol',
           pval_thresh=0.05,
           log2fc_thresh=0.75,
           fontsize=10, fontweight='bold',
           colors=['dimgrey', 'lightgrey', 'red'],
           top_right_frame=False,
           figsize=(5, 5), legend_pos=(1.4, 1),
           point_size=50,  # Added point_size parameter for fixed size
           save=False, 
           shapes=None,
           shape_order=None):
    
    '''
    Make a volcano plot from a pandas dataframe, with labeling for significant points and custom font options.
    
    data : pandas.DataFrame or path to csv
    log2fc : string, column name of log2 Fold-Change values
    pvalue : string, column name of the p-values to be converted to -log10 P values
    symbol : string, column name of gene IDs to use
    pval_thresh : numeric, threshold p-value for significance
    log2fc_thresh : numeric, threshold for log2 fold-change
    fontsize : int, size of labels
    fontweight : str, weight of font for labels (e.g., 'bold', 'normal')
    colors : list, colors to use (default: ['dimgrey', 'lightgrey', 'red'])
    top_right_frame : Boolean, Show the top and right frame (default: False)
    figsize : tuple, Size of figure (default: (5, 5))
    point_size : int, Fixed size for points (default: 50)
    save : boolean or string, option to save the plot
    '''
    
    if isinstance(data, str):
        df = pd.read_csv(data)
    else: 
        df = data.copy(deep=True)
        
    # Clean and imput 0s for p-value
    df = df.dropna()
    if df[pvalue].min() == 0:
        print('0s encountered for p value, imputing 1e-323')
        df[pvalue][df[pvalue] == 0] = 1e-323
        
    # Convert p-value threshold to -log10 scale and create a new column for plotting
    pval_thresh = -np.log10(pval_thresh)
    df['nlog10'] = -np.log10(df[pvalue])
    
    # Color light grey if below threshold, color DE red
    def map_color(a):
        log2FoldChange, nlog10 = a
        if abs(log2FoldChange) < log2fc_thresh or nlog10 < pval_thresh:
            return 'not DE'
        return 'DE'
    
    # Apply the map_color function
    df['color'] = df[[log2fc, 'nlog10']].apply(map_color, axis=1)
    
    # Set the order of colors: first 'DE' (for significant regulons), then 'not DE'
    hues = ['DE', 'not DE']
    
    # Set plot colors (DE -> red, not DE -> light grey)
    colors = colors[:len(df.color.unique())]
    
    # Create the volcano plot with fixed point size
    plt.figure(figsize=figsize)
    ax = sns.scatterplot(data=df, x=log2fc, y='nlog10',                 
                         hue='color', hue_order=hues, palette=colors,
                         s=point_size)  # Fixed point size for all points
    
    # Plot vertical and horizontal lines for thresholds
    ax.axhline(pval_thresh, zorder=0, c='k', lw=2, ls='--')
    ax.axvline(log2fc_thresh, zorder=0, c='k', lw=2, ls='--')
    ax.axvline(-log2fc_thresh, zorder=0, c='k', lw=2, ls='--')
    
    # Label significant regulons
    sig_df = df[(df['nlog10'] >= pval_thresh) & (abs(df[log2fc]) >= log2fc_thresh)]
    texts = []
    
    for i in range(len(sig_df)):
        txt = plt.text(sig_df.iloc[i][log2fc], 
                       sig_df.iloc[i]['nlog10'], 
                       sig_df.iloc[i][symbol], 
                       fontsize=fontsize, fontweight=fontweight)
        txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])
        texts.append(txt)
    
    # Adjust text to prevent overlap
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k', zorder=5))
    
    # Customize appearance
    for axis in ['bottom', 'left', 'top', 'right']:
        ax.spines[axis].set_linewidth(2)
        
    if not top_right_frame:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    ax.tick_params(width=2)
    plt.xticks(size=11, weight='bold')
    plt.yticks(size=11, weight='bold')
    plt.xlabel("$log_{2}$ fold change", size=15)
    plt.ylabel("-$log_{10}$ p-value", size=15)
    
    plt.legend(loc=1, bbox_to_anchor=legend_pos, frameon=False, prop={'weight':'bold'})
    
    # Save plot if required
    if save == True:
        files = os.listdir()
        for x in range(100):
            file_pref = "volcano_" + "%02d" % (x,)
            if len([x for x in files if x.startswith(file_pref)]) == 0:
                plt.savefig(file_pref + '.png', dpi=300, bbox_inches='tight')
                plt.savefig(file_pref + '.svg', bbox_inches='tight')
                break
    elif isinstance(save, str):
        plt.savefig(save + '.png', dpi=300, bbox_inches='tight')
        plt.savefig(save + '.svg', bbox_inches='tight')
                
    plt.show()

