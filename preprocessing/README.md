# Workflow of data preprocessing

1. Download data
   - For PBMC data provided by 10x, run `bash download_10x_PBMC.sh`
   - Other datasets are downloaded from links/GEO provided the authors. See Data availability statement in scMultiMap manuscript.
  
2. Call peaks by cell type
   - PBMC data: `preprocess_10x_PBMC_data.R`
   - brain data: `preprocess_brain_CG.R` and `preprocess_brain_ROSMAP.R`
  
3. MISC preprocessing analysis
   1. Merge 10x datasets on PBMC for more powerful analysis: `merge_10x_data.R`
   2. Save cell-type-specific objects for faster analysis: `save_seurat_obj_by_ct.R`
   3. Generate motif objects for trio analysis: `generate_motif_object.R`
   4. Generate candidate peak gene pairs between e.g. top 2000 highly expressed genes and top 20000 highly accessible peaks: `get_cis_peak_gene_pairs.R`. Replaced by `get_top_peak_gene_pairs` in scMultiMap R package.