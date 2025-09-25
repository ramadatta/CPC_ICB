# Regularly Requested Tasks 

## 1. Dotplot for a single gene across celltypes and conditions 


## 2. Check the signature across the conditions using specific set of genes using violin plots

```
    # 1. method1:
    sc.pl.violin(adata_Str, ["ACTG2", "ACTA2", "PLN"], groupby="ann_level_3_transferred_label")
    # 2. method2:  
    Analysis/1_Schiller_Lab/Projects/3_hPCLS_main_Combined_exvivo/hpcls_Exvivo_batch1_2.ipynb
        
    ### `Define Control and FC samples`
    
    `control_samples = ['D_zero_batch2', '12h_CC_batch2', 'D4_CC_batch1', 'D6_CC_batch1', 'D6_CC_batch2',] # I excluded the D_zero_batch1 since it was foundd to fibrotic genes
    fc_samples = ['12h_FC_batch1', '24h_FC_batch1','48h_FC_batch1','D4_FC_batch1','D6_FC_batch1',
    '12h_FC_batch2', '24h_FC_batch2','48h_FC_batch2','D4_FC_batch2','D6_FC_batch2']`
    
    ### `Create a new column in adata.obs to distinguish between Control and FC`
    
    `adata_Str.obs['group'] = adata_Str.obs['sample_batch'].apply(lambda x: 'Control' if x in control_samples else 'FC')`
    
    ### `Genes list for cell type "Fibroblasts"`
    
    `genes_of_interest = ['TCF21', 'PLIN2', 'CXCL1', 'ODC1', 'FGF10', 'COLEC12', 'EIF1', 'NFIL3']`
    
    `adata_subset = adata_Str[:, genes_of_interest]`
    
    ### `Group the data by the new 'group' column and calculate mean expression`
    
    `mean_expression_Str = adata_subset.to_df().groupby(adata_Str.obs['group']).mean()`
    
    # `Print the mean expression values`
    
    `print("Mean expression values:")
    print(mean_expression_Str)`
    
    # `Assume 'mean_expression_Str' is your DataFrame`
    
    `mean_expression_Str.reset_index(inplace=True)  # Ensure 'group' is a column and not an index
    long_format = pd.melt(mean_expression_Str, id_vars='group', var_name='Gene', value_name='Expression')`
    
    `print(long_format.head())`
    
    # `Assuming 'mean_expression_Str' is already calculated as shown above`
    
    # `Reset index to use 'group' as a column`
    
    `mean_expression_Str.reset_index(inplace=True)`
    
    # `Melt the DataFrame`
    
    `melted_data = pd.melt(mean_expression_Str, id_vars='group', var_name='Gene', value_name='Expression')`
    
    # `Calculate the overall mean expression for each group`
    
    `group_means = melted_data.groupby('group')['Expression'].mean().reset_index()`
    
    # `Create a violin plot`
    
    `plt.figure(figsize=(10, 6))
    sns.violinplot(x='group', y='Expression', data=melted_data)
    plt.title('Violin Plot of Mean Gene Expression Across Groups for SMC')
    plt.xlabel('Group')
    plt.ylabel('Mean Gene Expression')
    plt.show()`
```

### 3. Standardize the score for the UMAP plots
```   
    `sc.pp.scale(adata, max_value=10)`
```    

### 4. Line plots of gene expression across the timepoints 1 plot 1 gene (Refer 21 for a single gene)
```   
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from adjustText import adjust_text

def plot_gene_expression_over_time_sepPlot(
    adata, genes_of_interest, cluster_key, timepoints,
    plot_size=(8, 5), angle=90, colormap="tab10", label_fontsize=12
):
    """
    Plots the mean expression of a single gene across clusters in separate plots.
    Each gene is plotted individually with its own line.

    Parameters:
    - adata: AnnData object
    - genes_of_interest: list of str, gene names to plot (each will get its own plot)
    - cluster_key: str, the key in adata.obs containing cluster labels
    - timepoints: list of str, desired cluster order
    - plot_size: tuple, size of each plot (default is (8, 5))
    - angle: int, rotation angle for x-axis labels (default is 90)
    - colormap: str, matplotlib colormap name for gene colors (default is "tab10")
    - label_fontsize: int, font size for gene labels (default is 12)

    Returns:
    - None, but displays separate plots for each gene.
    """

    # Check for missing genes in adata
    missing_genes = [gene for gene in genes_of_interest if gene not in adata.var_names]
    if missing_genes:
        print(f"Warning: Missing genes in adata: {missing_genes}")

    # Filter for available genes
    valid_genes = [gene for gene in genes_of_interest if gene in adata.var_names]

    # Extract expression data and cluster labels
    expression_data = adata[:, valid_genes].to_df()
    expression_data["Cluster"] = adata.obs[cluster_key]

    # Calculate mean expression per cluster
    mean_expression = expression_data.groupby("Cluster").mean()

    # Filter and separate `CC` and `FC` clusters based on timepoints
    cc_clusters = [cluster for cluster in timepoints if "CC" in cluster and cluster in mean_expression.index]
    fc_clusters = [cluster for cluster in timepoints if "FC" in cluster and cluster in mean_expression.index]

    # Combine the ordered clusters: CC first, then FC
    ordered_clusters = cc_clusters + fc_clusters
    mean_expression = mean_expression.loc[ordered_clusters]

    # Map cluster labels to numeric positions for proper alignment
    cluster_positions = {cluster: idx for idx, cluster in enumerate(ordered_clusters)}

    # Generate unique colors for each gene using a colormap
    cmap = cm.get_cmap(colormap, len(valid_genes))
    gene_colors = {gene: cmap(i) for i, gene in enumerate(valid_genes)}

    # Loop through each gene and plot it separately
    for gene in valid_genes:
        plt.figure(figsize=plot_size)

        x_coords_cc = [cluster_positions[cluster] for cluster in cc_clusters]
        x_coords_fc = [cluster_positions[cluster] for cluster in fc_clusters]

        # Plot CC clusters with solid lines
        plt.plot(x_coords_cc, mean_expression.loc[cc_clusters, gene], marker="o", linestyle="-", color=gene_colors[gene], label="CC")

        # Plot FC clusters with dashed lines
        plt.plot(x_coords_fc, mean_expression.loc[fc_clusters, gene], marker="o", linestyle="--", color=gene_colors[gene], label="FC")

        # Add text annotations at the end of each line
        if cc_clusters:
            plt.text(x_coords_cc[-1] + 0.1, mean_expression.loc[cc_clusters[-1], gene],
                     f"{gene}", color=gene_colors[gene], fontsize=label_fontsize, va="center")

        if fc_clusters:
            plt.text(x_coords_fc[-1] + 0.1, mean_expression.loc[fc_clusters[-1], gene],
                     f"{gene}", color=gene_colors[gene], fontsize=label_fontsize, va="center")

        # Customize the plot
        plt.title(f"Mean Expression of {gene} Across Timepoints")
        plt.xlabel("Clusters")
        plt.ylabel("Mean Expression")
        plt.xticks(ticks=list(cluster_positions.values()), labels=ordered_clusters, rotation=angle)
        plt.legend()
        plt.tight_layout()

        # Show the plot
        plt.show()

``` 

### 5. Line plots of gene expression by each donor across the timepoints
```   
    Refer 21 for single gene
``` 

### 6. Matrix plot with mean z score
```
Other useful option is to normalize the gene expression using `sc.pp.scale`. Here we store this information under the `scale` layer. Afterwards we adjust the plot min and max and use a diverging color map (in this case `RdBu_r` where `_r` means reversed).

`# scale and store results in layer
pbmc.layers["scaled"] = sc.pp.scale(pbmc, copy=True).X`

`sc.pl.matrixplot(
    pbmc,
    marker_genes_dict,
    "clusters",
    dendrogram=True,
    colorbar_title="mean z-score",
    layer="scaled",
    vmin=-2,
    vmax=2,
    cmap="RdBu_r",
)`
```

### Note: Export gene signatures - **Analysis/1_Schiller_Lab/Projects/3_hPCLS_main_invivo/export_geneSignatures.ipynb**

### 7. Plot Gene Signature
```
genes_of_interest = ["TP63","CDH2"]

sc.tl.score_genes(adata, gene_list = genes_of_interest, score_name = "Basaloid_Score")
sc.pl.umap(adata, color = "Basaloid_Score", cmap = "Reds", size = 10)


sc.set_figure_params(figsize=(15, 10))
sc.pl.violin(adata, keys='Basaloid_Score', groupby='epi_harmony_leiden_cycle1_res200', rotation=90)
```

## 8. Check adata.var in a dictionary and assign a label if present in the dictionary
``` 
Convert the 'TF' column from human_tf into a set for quick lookup
tf_set = set(human_tf['TF'])
adata_invivo_Endo.var['Regulator'] = adata_invivo_Endo.var_names.map(lambda x: 'TF' if x in tf_set else None)
```

## 9. Relabel and create a new column
```
adata_exvivo_Endo_tf.obs['merged_controls'] = adata_exvivo_Endo_tf.obs['sample'].apply(lambda sample: 'Control' if sample == 'D_zero' or sample.endswith('_CC') else sample)

adata_exvivo_Epi.obs['fine_annot2'] = adata_exvivo_Epi.obs['fine_annot'].replace({'Aberrant Basaloid': 'Alveolar_epithelium', 
                                                                                                                 'AT1': 'Alveolar_epithelium', 
                                                                                                                 'AT2': 'Alveolar_epithelium'})
```

## 10. Create column by splitting
```
adata_exvivo_Endo_tf.obs[['Timepoint', 'Treated']] = adata_exvivo_Endo_tf.obs["merged_controls"].str.split('_',expand=True)
```
### 11. Filter gene list which are present in adata.var
```
filtered_upset_up_list = [gene for gene in upset_up_list if gene in adata_exvivo_Epi_tf.var_names]
```

### 12. Filter rows with NaN in pandas dataframe
```
filtered_df = df[df['MILO_IPFstages'].notnull()]
```
### 13. Exclude celltype(s)
```
adata_exvivo_Epi = adata_exvivo_Epi[~adata_exvivo_Epi.obs['fine_annot'].isin(['PNEC'])]
```
### 14. Just display particular celltype and rest of celltypes to be grey
```
# List of categories
categories = ["12h_CC", "12h_FC", "24h_FC", "48h_FC", "D4_FC", "D6_CC", "D6_FC", "D_zero"]

# Loop through categories and create UMAP plots
for category in categories:
    sc.pl.umap(
        adata_EpiCells,
        color=["sample"],
        groups=category,
        title=f'{category}',
        size=10
    )

```
### 15. extract DEGs from rank_gene_groups into dataframe
```
dedf = sc.get.rank_genes_groups_df(adata_Epi_Alv, group=None)
```
### 16. replace a portion of a string from a list 
```
random_regulons = [reg.replace('(+)','') for reg in random_regulons]
```
### 17. Note
```
In some cases, adata_Fibs.X might be a sparse matrix (e.g., scipy.sparse.csr_matrix). If that’s the case, you can convert it to a dense matrix before creating the DataFrame: adata_Fibs.X.toarray()
If adata_Fibs.X is already a dense numpy.ndarray, then there’s no need to convert it using adata_Fibs.X.toarray()
```

### 18. Plot Proportions
```
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Extract the color mapping from adata.uns
fine_annot2b_categories = adata.obs["fine_annot2b"].cat.categories
fine_annot2b_colors = adata.uns["fine_annot2b_colors"]

# Calculate proportions
condition_order = ['12h_CC','D4_CC', 'D6_CC', '12h_FC', '24h_FC', '48h_FC', 'D4_FC', 'D6_FC']
proportions = adata.obs.groupby(["sample", "fine_annot2b"]).size().unstack(fill_value=0)
proportions = proportions.div(proportions.sum(axis=1), axis=0)
proportions = proportions.reindex(condition_order, axis=0).fillna(0)

# Match colors to categories
colors = [fine_annot2b_colors[fine_annot2b_categories.get_loc(cat)] for cat in proportions.columns]

# Plot the stacked bar chart
plt.figure(figsize=(12, 8))
ax = proportions.plot(kind='bar', stacked=True, color=colors, ax=plt.gca())
plt.xlabel("Condition")
plt.ylabel("Proportion")
plt.title("Proportion of fine_annot2b across Condition Timepoints")
plt.xticks(rotation=90)

# Add proportion labels
for i, condition in enumerate(condition_order):
    cumulative_height = 0
    for j, category in enumerate(proportions.columns):
        value = proportions.loc[condition, category]
        if value > 0:
            ax.text(i, cumulative_height + value / 2, f'{value:.2f}', ha='center', va='center', fontsize=8)
        cumulative_height += value

# Customize the legend
plt.legend(
    title="fine_annot2b",
    labels=proportions.columns,
    loc='center left',
    bbox_to_anchor=(1.0, 0.5)
)
plt.grid(False)
plt.tight_layout()
plt.show()

```
### 19. Make gif from umaps
```
import scanpy as sc
import matplotlib.pyplot as plt
from PIL import Image
import os
import re

# List of categories
categories = ["Airway", "AT1", "AT2", "AT2 (SFTPBHi / CSF3RHi)", "AT2/Basaloid Transitional (TP63+/SFTPC+)", "Basaloid (CDH2+, TP63+)"]

# Create a directory to save the images
output_dir = "./umap_frames"
os.makedirs(output_dir, exist_ok=True)

# Function to sanitize file names
def sanitize_filename(name):
    return re.sub(r'[^\w\-_\. ]', '_', name)  # Replace non-alphanumeric characters with underscores

# Generate UMAP plots and save each as an image
image_files = []
for category in categories:
    sanitized_category = sanitize_filename(category)  # Sanitize the category name for the filename
    
    # Plot the UMAP
    sc.pl.umap(
        adata_Epi,
        color=["fine_annot1a"],
        groups=category,
        title=f'{category}',
        size=80,
        show=False,
        save=f"_{sanitized_category}.png"  # Save each plot with sanitized filename
    )
    
    # Move the saved plot to the output directory
    file_name = f"umap_{sanitized_category}.png"
    os.rename(f"figures/umap_{sanitized_category}.png", os.path.join(output_dir, file_name))
    image_files.append(os.path.join(output_dir, file_name))

# Create a GIF from the saved images
gif_output_path = "umap_plots_fine_annot1a.gif"
with Image.open(image_files[0]) as img:
    img.save(
        gif_output_path,
        save_all=True,
        append_images=[Image.open(file) for file in image_files[1:]],
        duration=500,  # Duration of each frame in milliseconds
        loop=0  # Loop infinitely
    )

print(f"GIF saved as {gif_output_path}")
```
### 20. Line plot from Gene Signature Score
```
import matplotlib.pyplot as plt
import pandas as pd

# Reorder the timepoints to ensure the correct order in the plot
ordered_timepoints = ['Control', '12h_FC', '24h_FC', '48h_FC', 'D4_FC', 'D6_FC']
data_df['Timepoint'] = pd.Categorical(data_df['Timepoint'], categories=ordered_timepoints, ordered=True)
data_df = data_df.sort_values('Timepoint')

# Calculate the mean score for each timepoint for both metrics
avg_scores = data_df.groupby('Timepoint').mean().reset_index()

# Plot both scores on the same line plot
plt.figure(figsize=(10, 6))

# Plotting Ectopic_EC_Score
plt.plot(avg_scores['Timepoint'], avg_scores['Ectopic_EC_Score'], marker='o', linestyle='-', color='blue', label='Ectopic EC Score')

# Plotting Ectopic_EC_2_3_DriverGeneScore
plt.plot(avg_scores['Timepoint'], avg_scores['Ectopic_EC_5_DriverGeneScore'], marker='o', linestyle='-', color='red', label='Ectopic EC 5 Driver Gene Score')

# Customize plot
plt.xlabel('Timepoint')
plt.ylabel('Average Score')
plt.title('Comparison of Ectopic EC Scores Across Timepoints')
plt.xticks(rotation=90)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
```
### 21. Line plot of Mean Gene Expression given a gene across multiple clusters
```
import pandas as pd
import matplotlib.pyplot as plt

def plot_gene_expression(adata, gene_of_interest, cluster_key, plot_size=(10, 6), line_color="blue"):
    """
    Plots the mean expression of a specified gene across clusters in a sorted line plot.

    Parameters:
    - adata: AnnData object
    - gene_of_interest: str, name of the gene to analyze
    - cluster_key: str, the key in adata.obs containing cluster labels
    - plot_size: tuple, size of the plot (default is (10, 6))
    - line_color: str, color of the line in the plot (default is "blue")
    
    Returns:
    - None, but displays a line plot
    """
    if gene_of_interest in adata.var_names:
        # Extract the cluster and gene expression data
        expression_data = adata[:, gene_of_interest].X
        cluster_data = adata.obs[cluster_key]

        # Create a DataFrame for easy manipulation
        df = pd.DataFrame({
            "Cluster": cluster_data,
            "Expression": expression_data.toarray().flatten() if hasattr(expression_data, "toarray") else expression_data.flatten()
        })

        # Calculate the mean expression per cluster
        mean_expression = df.groupby("Cluster")["Expression"].mean()

        # Sort clusters by mean expression
        mean_expression_sorted = mean_expression.sort_values()

        # Plot the sorted mean expression as a line plot
        plt.figure(figsize=plot_size)
        mean_expression_sorted.plot(kind="line", marker="o", color=line_color)
        plt.xticks(range(len(mean_expression_sorted)), mean_expression_sorted.index, rotation=0)
        plt.title(f"Mean Expression of {gene_of_interest} Across Clusters (Sorted by Expression)")
        plt.xlabel("Cluster (Sorted by Mean Expression, Lowest to Highest)")
        plt.ylabel("Mean Expression")
        plt.grid(False)
        plt.tight_layout()  # Adjust layout for better fit
        plt.show()
    else:
        print(f"Gene {gene_of_interest} not found in adata.")

plot_gene_expression(adata, gene_of_interest="TP63", cluster_key="epi_harmony_leiden_cycle1_res200", plot_size=(10, 6), line_color="green")
```
### 22. Line plot of Mean Gene Expression given set of gene across multiple timepoints in order
```
def plot_combined_gene_lines_custom_order_adjusted_labels(adata, genes_of_interest, cluster_key, timepoints,
                                                          plot_size=(12, 6), cc_color="blue", fc_color="red",
                                                         angle=90):
    """
    Plots the mean expression of multiple genes across clusters in a line plot,
    with `CC` clusters appearing first on the x-axis, followed by `FC`, in a custom order.
    Annotates gene names at the end of the lines with adjusted positions using `adjustText`.
    
    Parameters:
    - adata: AnnData object
    - genes_of_interest: list of str, gene names to plot
    - cluster_key: str, the key in adata.obs containing cluster labels
    - timepoints: list of str, desired cluster order
    - plot_size: tuple, size of the plot (default is (12, 6))
    - cc_color: str, color for CC groups (default is "blue")
    - fc_color: str, color for FC groups (default is "red")
    
    Returns:
    - None, but displays a line plot
    """
    # Check for genes in adata
    missing_genes = [gene for gene in genes_of_interest if gene not in adata.var_names]
    if missing_genes:
        print(f"Warning: Missing genes in adata: {missing_genes}")

    # Filter for available genes
    valid_genes = [gene for gene in genes_of_interest if gene in adata.var_names]

    # Extract expression data and cluster labels
    expression_data = adata[:, valid_genes].to_df()
    expression_data["Cluster"] = adata.obs[cluster_key]

    # Calculate mean expression per cluster
    mean_expression = expression_data.groupby("Cluster").mean()

    # Filter and separate `CC` and `FC` clusters based on timepoints
    cc_clusters = [cluster for cluster in timepoints if "CC" in cluster and cluster in mean_expression.index]
    fc_clusters = [cluster for cluster in timepoints if "FC" in cluster and cluster in mean_expression.index]

    # Combine the ordered clusters: CC first, then FC
    ordered_clusters = cc_clusters + fc_clusters
    mean_expression = mean_expression.loc[ordered_clusters]

    # Map cluster labels to numeric positions for proper alignment
    cluster_positions = {cluster: idx for idx, cluster in enumerate(ordered_clusters)}

    # Initialize the plot
    plt.figure(figsize=plot_size)

    # Collect texts for adjustment
    texts = []

    # Plot CC clusters
    for gene in valid_genes:
        x_coords = [cluster_positions[cluster] for cluster in cc_clusters]
        plt.plot(x_coords, mean_expression.loc[cc_clusters, gene], marker="o", linestyle="-", color=cc_color)
        # Add text annotations for CC clusters
        texts.append(plt.text(x_coords[-1] + 0.1, mean_expression.loc[cc_clusters[-1], gene], f"{gene}", color=cc_color, fontsize=10, va="center"))

    # Plot FC clusters
    for gene in valid_genes:
        x_coords = [cluster_positions[cluster] for cluster in fc_clusters]
        plt.plot(x_coords, mean_expression.loc[fc_clusters, gene], marker="o", linestyle="--", color=fc_color)
        # Add text annotations for FC clusters
        texts.append(plt.text(x_coords[-1] + 0.1, mean_expression.loc[fc_clusters[-1], gene], f"{gene}", color=fc_color, fontsize=10, va="center"))

    # Adjust text to prevent overlap
    adjust_text(texts, arrowprops=dict(arrowstyle="->", color="gray", lw=0.5))

    # Customize the plot
    plt.title("Mean Expression of Genes Across Timepoints")
    plt.xlabel("Clusters")
    plt.ylabel("Mean Expression")
    plt.xticks(ticks=list(cluster_positions.values()), labels=ordered_clusters, rotation=angle)
    plt.tight_layout()
    plt.show()


plot_combined_gene_lines_custom_order_adjusted_labels(
    adata_Endo,
    genes_of_interest=gene_of_int_CC,
    cluster_key="sample",
    timepoints=timepoints
)
```
### 23. Violin plot for subset clusters
```
# Specify the two cell types to compare
celltypes_to_compare = ["Ectopic EC (PLVAP/VWA1)", "Ectopic EC (VWF)"]

# Subset the AnnData object to include only these cell types
adata_temp = adata_Endo[adata_Endo.obs["fine_annot4b"].isin(celltypes_to_compare)]

# Plot the violin plot
sc.pl.violin(
    adata_temp,
    gene_of_int_CC, ## Ectopic EC in CC 
    groupby="fine_annot4b"
)
```
### 24. Melt a matrix (Wide) to Long Format
```
# Reshape the wide DataFrame into long format
ulm_estimate_long = adata.obsm['ulm_estimate'].reset_index().melt(
    id_vars='index', var_name='Gene', value_name='Value'
)

# Rename the `index` column to `Cell`
ulm_estimate_long.rename(columns={'index': 'Cell'}, inplace=True)
```

### 25. Check the sizes of adata components and layers
```
from scipy.sparse import issparse

# Check if adata.X is sparse
if issparse(adata_wo_Cil_Deu.X):
    X_nbytes = (
        adata_wo_Cil_Deu.X.data.nbytes +
        adata_wo_Cil_Deu.X.indices.nbytes +
        adata_wo_Cil_Deu.X.indptr.nbytes
    )
    print(f"X: {X_nbytes / 1e9:.2f} GB (sparse matrix)")
else:
    print(f"X: {adata_wo_Cil_Deu.X.nbytes / 1e9:.2f} GB (dense matrix)")

# Check the size of obs
obs_nbytes = adata_wo_Cil_Deu.obs.memory_usage(deep=True).sum()
print(f"obs: {obs_nbytes / 1e6:.2f} MB")

# Check the size of var
var_nbytes = adata_wo_Cil_Deu.var.memory_usage(deep=True).sum()
print(f"var: {var_nbytes / 1e6:.2f} MB")

# Check the size of uns
import sys

uns_nbytes = sum(sys.getsizeof(adata_wo_Cil_Deu.uns[key]) for key in adata_wo_Cil_Deu.uns.keys())
print(f"uns: {uns_nbytes / 1e6:.2f} MB")

# Check the size of layers
layers_nbytes = 0
for layer_name, layer_data in adata_wo_Cil_Deu.layers.items():
    if issparse(layer_data):
        layer_nbytes = (
            layer_data.data.nbytes +
            layer_data.indices.nbytes +
            layer_data.indptr.nbytes
        )
    else:
        layer_nbytes = layer_data.nbytes
    layers_nbytes += layer_nbytes
print(f"layers: {layers_nbytes / 1e9:.2f} GB")


# Check the size of each layer in adata
from scipy.sparse import issparse

# Iterate over layers to compute their sizes
for layer_name, layer in adata_wo_Cil_Deu.layers.items():
    if issparse(layer):
        print(f"Layer {layer_name}: {layer.data.nbytes / 1e9:.2f} GB (sparse)")
    else:
        print(f"Layer {layer_name}: {layer.nbytes / 1e9:.2f} GB (dense)")
```

### 26. Seurat Style Ranking of Markers from scanpy

```
import numpy as np
import pandas as pd

def get_markers(
    adata,
    groupby,
    key="rank_genes_groups",
    p_val_cutoff=0.05,
    logfc_cutoff=0.5,
    pct_1_cutoff=0.5,
    pct_2_cutoff=None,
    ranking_metric=True,
    layer="counts"  # Layer to specify the data source
):
    """\
    Extract markers using a specified adata layer for raw counts with gene name handling.
    """
    # Extract markers and filter by logFC and p-value cutoffs
    markers = pd.concat([
        pd.DataFrame(adata.uns[key]["names"]).melt(),
        pd.DataFrame(adata.uns[key]["pvals_adj"]).melt(),
        pd.DataFrame(adata.uns[key]["logfoldchanges"]).melt()
    ], axis=1)
    markers.columns = ("cluster", "gene", "cluster2", "p_val_adj", "cluster3", "avg_logFC")
    markers = markers.loc[:, ["cluster", "gene", "avg_logFC", "p_val_adj"]]
    markers = markers.loc[markers.avg_logFC > logfc_cutoff, ]
    markers = markers.loc[markers.p_val_adj < p_val_cutoff, ]
    markers["pct.1"] = pd.Series(dtype=float)
    markers["pct.2"] = pd.Series(dtype=float)

    # Select the data matrix (layer or main matrix)
    data_matrix = adata.layers[layer] if layer else adata.X

    # Calculate pct.1 and pct.2 for each gene in each cluster
    for cluster in markers.cluster.unique():
        cells = adata.obs[groupby] == cluster
        in_cluster_selector = markers.cluster == cluster
        genes = markers.gene[in_cluster_selector]
        
        # Convert gene names to indices
        gene_indices = adata.var_names.get_indexer(genes)
        
        # Calculate pct.1 using gene indices
        in_cluster = np.sum(data_matrix[cells, :][:, gene_indices].toarray() > 0, axis=0) / cells.sum()
        markers.loc[in_cluster_selector, "pct.1"] = in_cluster.T

        # Calculate pct.2 using gene indices
        other_cells = adata.obs[groupby] != cluster
        other_clusters = np.sum(data_matrix[other_cells, :][:, gene_indices].toarray() > 0, axis=0) / other_cells.sum()
        markers.loc[in_cluster_selector, "pct.2"] = other_clusters.T

    # Filter by pct.1 and pct.2 cutoffs
    markers = markers.loc[markers["pct.1"] > pct_1_cutoff, ]
    if pct_2_cutoff is not None:
        markers = markers.loc[markers["pct.2"] < pct_2_cutoff, ]

    # Calculate ranking score if requested and data exists
    if ranking_metric and not markers.empty and "pct.1" in markers and "pct.2" in markers:
        markers["score"] = markers["avg_logFC"] * (markers["pct.1"] / (markers["pct.2"] + 1e-10))
        markers = markers.sort_values(by="score", ascending=False)

    markers["p_val"] = markers.p_val_adj
    # Ensure 'score' column exists for consistency
    if "score" not in markers.columns:
        markers["score"] = np.nan
    markers = markers.loc[:, ["p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene", "score"]]
    
    return markers
```
### 27. Annotate clusters in scanpy
```
types = {
    "AT1": ["16","26", "17", "11", "1","21"], # 17 with pinch of salt
     "AT2-AT1 differentiating": ["26"],
    "Basaloid Transitional":[ "5", "9"],
    "AT2 (TRB-SC/AT0-like)" : ["20","4","12"],
    "Basaloid": ["3", "19"],
    "AT2": ["0", "2", "6","14", "15","25","7", "10", "18", "22", "24"],  
    "Airway" : ["8", "13", "23"],
    # "Club": ["3", "4"],
    # "Goblet": ["3", "4"],
    # "Ciliated": ["23"],
}


# Create a dictionary to map cluster numbers to cell types
cluster_to_cell_type = {cluster_num: cell_type for cell_type, cluster_list in types.items() for cluster_num in cluster_list}

# Create a new column in adata.obs with cell types based on leiden_res2
adata.obs["fine_annot2a"] = [cluster_to_cell_type[i] for i in adata.obs["epi_harmony_leiden_cycle1_res200"]]
adata.obs["fine_annot2a"].value_counts()
```

### 28. Matrixplot Rank Gene Groups
```
sc.tl.rank_genes_groups(adata, "fine_annot2a", method="wilcoxon")
sc.pl.rank_genes_groups_matrixplot(
    adata_Epi_Scenic,
    n_genes=5,
    key="Scenic_Regulon_markers_fine_annot2b",
    groupby = "fine_annot2b",
    standard_scale="var"
)

```

### 29. Top genes from scanpy dataframe
```
top10_by_celltype = drdf.groupby("group").apply(lambda x: x.nlargest(10, "scores")).reset_index(drop=True)
```

### 30. Save high quality images from Scanpy
```
# Overall celltype
sc.set_figure_params(figsize=(10,10), dpi=150, dpi_save=600)
# This one
markers = ['NKX2-1']
sc.pl.umap(adata_siRNA_epi, color=markers, cmap=grey_red, save="NKX2-1.png", size=60)
```
### 31. Extract gene signatures from Scanpy
```
dedf = sc.get.rank_genes_groups_df(adata, group=None)
top_100_per_group = dedf.groupby("group").head(100)[["group", "names"]]
top_100_per_group.to_csv("invivo_immune_top_100_per_group.csv", index=False)
```

### 32. Subcluster then run this
```
%%time
print("HVG")
sc.pp.highly_variable_genes(adata_AT1)
print("PCA")
sc.pp.pca(adata_AT1, n_comps=50, use_highly_variable=True, svd_solver="arpack")
print("Neighbors")
sc.pp.neighbors(adata_AT1, use_rep="X_pca_harmony")
print("UMAP")
sc.tl.umap(adata_AT1)
print("leiden")
sc.tl.leiden(adata_AT1, resolution = 0.5, key_added = "AT1_leiden_res50")
```

### 33. Bigger R plots in JupyterLab
```
options(repr.plot.width = 15, repr.plot.height = 10, repr.plot.res = 300)
```

### 34. Automatically check the genes present in adata or not and plot only those exist
```
var_names=[g for g in genelist if g in adata.var_names],
```

### 35. Seurat object to h5ad conversion
```
# This is an alternative to sceasy
library(Seurat)
library(SeuratDisk)
seurat_obj <- readRDS("GSE135893_ILD_annotated_fullsize.rds")
SaveH5Seurat(seurat_obj, filename = "GSE135893_ILD_annotated_fullsize.h5Seurat")
seurat_obj <- UpdateSeuratObject(seurat_obj)
SaveH5Seurat(seurat_obj, filename = "GSE135893_ILD_annotated_fullsize.h5Seurat")
Convert("GSE135893_ILD_annotated_fullsize.h5Seurat", dest = "h5ad", overwrite = TRUE)
```
