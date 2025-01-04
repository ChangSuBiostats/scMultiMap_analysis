# From https://rdrr.io/cran/Seurat/src/R/preprocessing.R

# Convert SCTModel class to vst_out used in the sctransform
# @param SCTModel
# @return Return a list containing sct model
#
SCTModel_to_vst <- function(SCTModel) {
  feature.params <- c("theta", "(Intercept)",  "log_umi")
  feature.attrs <- c("residual_mean", "residual_variance" )
  vst_out <- list()
  vst_out$model_str <- slot(object = SCTModel, name = "model")
  vst_out$model_pars_fit <- as.matrix(x = slot(object = SCTModel, name = "feature.attributes")[, feature.params])
  vst_out$gene_attr <- slot(object = SCTModel, name = "feature.attributes")[, feature.attrs]
  vst_out$cell_attr <- slot(object = SCTModel, name = "cell.attributes")
  vst_out$arguments <- slot(object = SCTModel, name = "arguments")
  return(vst_out)
}


