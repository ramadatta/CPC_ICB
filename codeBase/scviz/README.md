# scviz

Single-cell visualization tools for creating interactive plots with UMAP and proportional analysis.

## Features

- Interactive UMAP visualization with cell type highlighting
- Stacked bar charts showing cell type proportions across conditions
- Toggle between proportions and absolute counts
- Support for single or multiple condition comparisons
- Exportable as standalone HTML files

## Functions

### `interactive_umap_propbar_plot()`

Creates an interactive plot combining UMAP visualization with stacked bar charts for cell type proportion analysis across different conditions.

Key capabilities:
- Visualize cell distributions in UMAP space
- Compare cell type proportions across conditions
- Highlight specific cell types across both UMAP and bar charts
- Switch between proportion and count views
- Support for one or two condition groupings (tabs)
- Export as shareable HTML files like [this](https://drive.google.com/file/d/15dA1TPUG5N2OgfECSGMqc0cZf7cR5wzV/view?usp=drive_link)
  
```
interactive_umap_propbar_plot(
    adata=adata_Epi,
    celltype_col='cytetype_annotation_epi_harmony_leiden_cycle1_res300',
    color_key='cytetype_annotation_epi_harmony_leiden_cycle1_res300_colors',
    condition_col_tab1='condition_1c_Timepoints_Control',
    condition_order_tab1=['Control', '12h_FC', '24h_FC', '48h_FC', 'D4_FC', 'D6_FC'],
    tab1_title='CC vs Timepoints',
    condition_col_tab2='condition_1b_FC_Control',
    condition_order_tab2=['Control', 'FC'],
    tab2_title='CC vs FC',
    output_file='testing_both_interactive_cell_proportions_with_tabs.html'
)

# Example usage with only tab1:
interactive_umap_propbar_plot(
    adata=adata_Epi,
    celltype_col='cytetype_annotation_epi_harmony_leiden_cycle1_res300',
    color_key='cytetype_annotation_epi_harmony_leiden_cycle1_res300_colors',
    condition_col_tab1='condition_1c_Timepoints_Control',
    condition_order_tab1=['Control', '12h_FC', '24h_FC', '48h_FC', 'D4_FC', 'D6_FC'],
    tab1_title='CC vs Timepoints',
    output_file='testing_tab1_interactive_cell_proportions_single_tab.html'
)
```
