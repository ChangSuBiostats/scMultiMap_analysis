# scMultiMap_analysis

This GitHub repository contains the code for reproducing the results in 

Chang Su, Dongsoo Lee, Peng Jin and Jingfei Zhang. (2024). [Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data](https://www.biorxiv.org/content/10.1101/2024.09.24.614814v1). Manuscript under revision at *Nature Communications*.

## Structure

- method_implementation: implementations of Signac and SCENT
- preprocessing: code for data processing 
- simulation: code for evaluating type I errors and power based on simulation & permutation
- real_data_analysis: code for analysis based on real data
- manuscript: code for generating the figures

Please read **workflow.md** or **README.md** under each folder for detailed documentation.

## R package and vignettes
Please visit https://github.com/ChangSuBiostats/scMultiMap for the R package that implements `scMultiMap`, and https://changsubiostats.github.io/scMultiMap/ for vignettes.
