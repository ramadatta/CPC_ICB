volcano(results, log2fc='log2FC', pvalue='pval', symbol='regulon', log2fc_thresh=0.5, 
        colors = ['red', 'blue', 'red'], pval_thresh=0.05, fontsize=8, fontweight='bold', point_size=10, figsize=(7, 7))
        
        
pSTY_gene_list = ["EIF4EBP1", "EP300", "FAM195B", "FLNA", "FOXO3", "MAF1", "MAPK1", "MAPK14", "MAPK3", 
              "NFATC4", "PKAR2B", "PXN", "RB1", "SSB", "SYMPK", "TP53", "VIM", "YBX1"]
              
# Example list of manually labeled genes
manual_gene_list = pSTY_gene_list

# Example gene color dictionary
# Assign a single color to all genes in the list (e.g., red)
gene_color_dict = {gene: "blue" for gene in pSTY_gene_list}

# Create the volcano plot, manually labeling the genes and coloring them
volcano_pass_genelist(
    results, 
    log2fc='log2FC', 
    pvalue='pval', 
    symbol='regulon',
    log2fc_thresh=0.5, 
    pval_thresh=0.05, 
    point_size=6,
    manual_labels=manual_gene_list,  # Manually label these genes
    gene_color_dict=gene_color_dict,  # Manually color these genes
    figsize=(6, 6), colors = ['grey', 'red', 'red'],
    fontsize=7, fontweight='bold'
)              
