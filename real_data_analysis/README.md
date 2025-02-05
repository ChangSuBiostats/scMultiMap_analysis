# Workflow of real_data_analyais

1. Evaluate scMultiMap, SCENT and Signac on PBMC data and brain data for the analysis in 2-5.
  - `bash analysis.sh`
2. Reproducibility analysis
  - reproducibility_PBMC.ipynb
  - reproducibility_brain.ipynb
3. Consistency with orthogonal modalities
  - validation_PBMC.ipynb
  - validation_brain_CG_plac_seq.ipynb
4. Trio analysis
  - trio_PBMC.ipynb
  - trio_brain.ipynb
5. Benchmarking running time
  - benchmark_computing_time.R
6. MISC analysis
  - validation_brain_CG_ENCODE.ipynb
  - LDSC analysis: folder LDSC/
7. Other analysis covered by scMultiMap R packages
  - Differential co-expression: Please refer to scMultiMap R package, vignette [scMultiMap for disease-control studies](https://changsubiostats.github.io/scMultiMap/articles/disease_control.html)
  - Integrative analysis with GWAS results: Please refer to scMultiMap R package, vignette [scMultiMap for integrative analysis with GWAS results](https://changsubiostats.github.io/scMultiMap/articles/GWAS.html)
8. Helper functions
  - reproducibility_helper.R
  - validation_helper.R
  - trio_anlaysis_helper.R
  - peak_annotation_helper.R
