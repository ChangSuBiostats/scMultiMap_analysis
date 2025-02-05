# README for functions that implements scMultiMap, SCENT and Signac 

1. `moment_regressions.R`
   - An archived implementation of scMultiMap used during method development
   - Gives exactly the same numerical results as implemented in R package
    
2. `SCENT_early_stop.R`
   - An implementation of SCENT that fixed the number of boostrap replicates to avoid extreme computational costs
   - Based on its original implementation at https://github.com/immunogenomics/SCENT/blob/main/R/SCENTfunctions.R

3. `Signac_LinkPeaks.R`
   - A customized implementation of Signac that fixes the coding error in p value calculation and speeds up computation
   - Based on its original implementation at https://github.com/stuart-lab/signac/blob/8ecdde2/R/links.R

For more details on the customization made in 2-3, please refer to the Methods section in the manuscript.
