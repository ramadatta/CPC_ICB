# Frequently asked tasks

## 1. Dotplot for a single gene across celltypes and conditions 


## 2. Check the signature across the conditions using specific set of genes using violin plots

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
    
3. Standardize the score for the UMAP plots
    
    `sc.pp.scale(adata, max_value=10)`
    
4. Line plots of gene expression acroos the timepoints
5. Line of gene expression by each donor across the timepoints
6. Matrix plot with mean z score

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

1. Export gene signatures - **Analysis/1_Schiller_Lab/Projects/3_hPCLS_main_invivo/export_geneSignatures.ipynb**
2. Plot Gene Signature

# `Genes list for cell type "SMC" (Smooth Muscle Cells)`

`genes_of_interest = adj_str_df_NFATC4_top20.target.tolist()`

`sc.tl.score_genes(adata_STR_wo_EMPH, gene_list = genes_of_interest, score_name = "NFATC4_TargetGene_Score")`

`sc.pl.umap(adata_STR_wo_EMPH, color = "NFATC4_TargetGene_Score", cmap = "Reds", size = 10)`

1. we can also draw violin plot using the gene signature

sc.set_figure_params(figsize=(15, 10))

# Assuming 'cell_type' is the column in adata.obs that categorizes cells into different types

# Ensure 'NFATC4_TargetGene_Score' is in adata.obs and 'cell_type' is correctly set up

sc.pl.violin(adata_STR_wo_EMPH, keys='NFATC4_TargetGene_Score', groupby='SEACell_predom_cellstate',
rotation=90)


## 3. Check adata.var in a dictionary and assign a label if present in the dictionary
``` 
Convert the 'TF' column from human_tf into a set for quick lookup
tf_set = set(human_tf['TF'])
adata_invivo_Endo.var['Regulator'] = adata_invivo_Endo.var_names.map(lambda x: 'TF' if x in tf_set else None)
```

## 4. Relabel and create a new column
```
adata_exvivo_Endo_tf.obs['merged_controls'] = adata_exvivo_Endo_tf.obs['sample'].apply(lambda sample: 'Control' if sample == 'D_zero' or sample.endswith('_CC') else sample)

adata_exvivo_Epi.obs['fine_annot2'] = adata_exvivo_Epi.obs['fine_annot'].replace({'Aberrant Basaloid': 'Alveolar_epithelium', 
                                                                                                                 'AT1': 'Alveolar_epithelium', 
                                                                                                                 'AT2': 'Alveolar_epithelium'})
```

## 5. Create column by splitting
```
adata_exvivo_Endo_tf.obs[['Timepoint', 'Treated']] = adata_exvivo_Endo_tf.obs["merged_controls"].str.split('_',expand=True)
```
### 6. Filter gene list which are present in adata.var
```
filtered_upset_up_list = [gene for gene in upset_up_list if gene in adata_exvivo_Epi_tf.var_names]
```

### 7. Filter rows with NaN in pandas dataframe
```
filtered_df = df[df['MILO_IPFstages'].notnull()]
```
### 8. Exclude celltype(s)
```
adata_exvivo_Epi = adata_exvivo_Epi[~adata_exvivo_Epi.obs['fine_annot'].isin(['PNEC'])]
```
### 9. Just display particular celltype and rest of celltypes to be grey
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
### 10. extract DEGs from rank_gene_groups into dataframe
```
dedf = sc.get.rank_genes_groups_df(adata_Epi_Alv, group=None)
```
### 11. replace a portion of a string from a list 
```
random_regulons = [reg.replace('(+)','') for reg in random_regulons]
```
### 12. Note
```
In some cases, adata_Fibs.X might be a sparse matrix (e.g., scipy.sparse.csr_matrix). If that’s the case, you can convert it to a dense matrix before creating the DataFrame: adata_Fibs.X.toarray()
If adata_Fibs.X is already a dense numpy.ndarray, then there’s no need to convert it using adata_Fibs.X.toarray()
```
