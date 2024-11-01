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

# TO RUN
# markers = get_markers(
#     adata,
#     groupby="leiden_cluster",
#     p_val_cutoff=0.1,
#     logfc_cutoff=0.5,
#     pct_1_cutoff=0.3,
#     pct_2_cutoff=0.1,          # Optional: less than 50% in other groups
#     ranking_metric=True,
#     layer="counts"             # Specify the layer to use for raw counts
# )

# p_val	avg_logFC	pct.1	pct.2	p_val_adj	cluster	gene	score
# 4	4.251024e-199	9.344220	0.660668	0.002902	4.251024e-199	Lymphatic EC	SEMA3D	2127.518516
# 35	1.727526e-91	8.459714	0.443445	0.002394	1.727526e-91	Lymphatic EC	SV2B	1567.068526
# 28	2.976919e-103	8.922918	0.467866	0.002974	2.976919e-103	Lymphatic EC	PROX1	1403.626643
# 7	1.319426e-172	9.541409	0.619537	0.006094	1.319426e-172	Lymphatic EC	FLRT2	970.079783
# 16	2.505697e-121	7.837757	0.510283	0.005151	2.505697e-121	Lymphatic EC	RELN	776.517231
# 93061	1.657343e-98	1.794000	0.301198	0.084486	1.657343e-98	Arterial EC	INPP5D	6.395748
# 93059	2.468627e-100	1.754441	0.306886	0.091209	2.468627e-100	Arterial EC	LXN	5.903081
