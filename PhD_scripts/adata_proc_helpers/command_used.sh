time python /home/ramadatta/sw/scripts/adata_proc_helpers/get_value_counts.py adata_healthy_IPF_fromHLCA_LisaSikkema_2_afterHarmony.h5ad study
Value counts for 'study':
Kaminski_2020             239707
Banovich_Kropski_2020     179562
Sheppard_2020              61603
Misharin_Budinger_2018     56737


python /home/ramadatta/sw/scripts/adata_proc_helpers/subset_adata.py adata_healthy_IPF_fromHLCA_LisaSikkema_2_afterHarmony.h5ad study Sheppard_2020 Sheppard_2020_healthy_IPF.h5ad

time python /home/ramadatta/sw/scripts/adata_proc_helpers/subset_adata.py GLPG_v022024_Stromal_annotated_cellstates_v1.h5ad annotation_level2b "Stromal - fibroblast" GLPG_v022024_Stromal_annotated_cellstates_v1_Fibs.h5ad &

python /home/ramadatta/sw/scripts/adata_proc_helpers/get_unique_values_from_obs.py adata_healthy_IPF_fromHLCA_LisaSikkema_2_afterHarmony.h5ad study

python /home/ramadatta/sw/scripts/adata_proc_helpers/show_adcols.py Kaminski_2020_healthy_IPF.h5ad
