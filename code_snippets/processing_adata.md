# Transfer labels from one adata to another

```
import pandas as pd
# Extract the relevant columns from adata_scarches
scarches_labels = adata_scarches.obs[[
    'ann_level_1_transferred_label', 
    'ann_level_2_transferred_label', 
    'ann_level_3_transferred_label', 
    'ann_level_4_transferred_label', 
    'ann_level_5_transferred_label'
]]

# Ensure the index of scarches_labels is the cell barcodes
scarches_labels.index = adata_scarches.obs_names

# Align the cell barcodes between adata and adata_scarches
common_barcodes = adata.obs_names.intersection(scarches_labels.index)

# Subset both adata and scarches_labels to the common barcodes
adata_subset = adata[common_barcodes, :]
scarches_labels_subset = scarches_labels.loc[common_barcodes, :]

# Transfer the labels from one adata to another adata
adata_subset.obs[[
    'ann_level_1_transferred_label', 
    'ann_level_2_transferred_label', 
    'ann_level_3_transferred_label', 
    'ann_level_4_transferred_label', 
    'ann_level_5_transferred_label'
]] = scarches_labels_subset

# If you want to update the original adata with the transferred labels
adata.obs.loc[common_barcodes, [
    'ann_level_1_transferred_label', 
    'ann_level_2_transferred_label', 
    'ann_level_3_transferred_label', 
    'ann_level_4_transferred_label', 
    'ann_level_5_transferred_label'
]] = scarches_labels_subset
```

# Exclude cells by counts
```
# Calculate the counts of each label in ann_level_3_transferred_label
label_counts = adata.obs['ann_level_3_transferred_label'].value_counts()

# Identify labels with counts >= 100
valid_labels = label_counts[label_counts >= 100].index

# Filter adata to exclude cells with labels having counts < 100
adata_filtered = adata[adata.obs['ann_level_3_transferred_label'].isin(valid_labels)]
```

# Exclude rows by pattern
```
# Exclude cells where ann_level_3_transferred_label is "Unknown"
adata_subset = adata_EPI_Scenic[adata_EPI_Scenic.obs['ann_level_3_transferred_label'] != 'Unknown']
```
# Replace a value in adata.obs column with another
```
adata.obs["batch"]=adata.obs["batch"].replace(['0'], 'Batch1')
```
# tabulate counts by two columns
```
pd.crosstab(adata.obs["batch"], adata.obs["sample"])
```
