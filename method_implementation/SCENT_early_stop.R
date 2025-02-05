# Reference: https://github.com/immunogenomics/SCENT/blob/main/R/SCENTfunctions.R
SCENT_algorithm_early_stop <- function(object, celltype, ncores, regr = 'poisson', bin = TRUE,
                                      early_stop_R = 700, fixed=F, seed=2023){
    
    set.seed(2023) ##-(Chang) for reproducibility
  res <- data.frame()
  for (n in 1:nrow(object@peak.info)){ ####c(1:nrow(chunkinfo))
    gene <- object@peak.info[n,1] #GENE is FIRST COLUMN OF PEAK.INFO
    this_peak <- object@peak.info[n,2] #PEAK is SECOND COLUMN OF PEAK.INFO
      print(n)
    atac_target <- data.frame(cell = colnames(object@atac), atac = object@atac[this_peak,])
    #binarize peaks:
    if(bin){
        if(sum(atac_target$atac>0) > 0) { ### new
      atac_target[atac_target$atac>0,]$atac<-1
            }
    }
    mrna_target <- object@rna[gene,]
    df <- data.frame(cell=names(mrna_target),exprs=as.numeric(mrna_target))
    df<-merge(df,atac_target,by='cell')
      
    df<-merge(df,object@meta.data,by='cell')
    df2 <- df[df[[object@celltypes]] == celltype,]
    nonzero_m  <- length( df2$exprs[ df2$exprs > 0] ) / length( df2$exprs )
    nonzero_a  <- length( df2$atac[ df2$atac > 0] ) / length( df2$atac )
      if(nonzero_a == 0){
          print('return NA for peaks==0')
          outs <- data.frame(gene = gene, peak = this_peak, beta=NA, se=NA, z=NA, p=NA,boot_basic_p=NA)
          res<-rbind(res,out)
          next
        }
      #print(nonzero_m)
      #print(nonzero_a)
    #if(nonzero_m > 0.05 & nonzero_a > 0.05){
      ##-(Chang) remove thresholding by sparsity, as we are already taking the top 2k and 20k genes and peaks
      if(T){
      #Run Regression Once Before Bootstrapping:
      res_var <- 'exprs'
      pred_var <- c('atac', object@covariates)
      formula <- as.formula(paste(res_var, paste(pred_var, collapse = '+'), sep = '~'))
      #Estimated Coefficients Obtained without Bootstrapping:
      if(regr == 'poisson'){
        base = glm(formula, family = 'poisson', data = df2)
        coefs<-summary(base)$coefficients['atac',]
        assoc <- assoc_poisson
      } else if (regr == 'negbin'){
        base = glm.nb(formula, data = df2)
        coefs<-summary(base)$coefficients['atac',]
        assoc <- assoc_negbin
      }
      ###Iterative Bootstrapping Procedure: Estimate the Beta coefficients and associate a 2-sided p-value.
      bs = boot::boot(df2,assoc, R = 100, formula = formula, stype = 'i', parallel = 'multicore', ncpus = ncores)
      p0 = basic_p(bs$t0[1], bs$t[,1])
      if(!fixed){
          #  old_bs_t <- bs$t[,1]
      if(p0<0.1){
        bs = boot::boot(df2,assoc, R = 100, formula = formula,  stype = 'i', parallel = 'multicore', ncpus = ncores)
        p0 = basic_p(bs$t0[1], bs$t[,1])
      #    old_bs_t <- c(old_bs_t, bs$t[,1])
      }
      if(p0<0.05){
        bs = boot::boot(df2,assoc, R = early_stop_R, formula = formula,  stype = 'i', parallel = 'multicore', ncpus = ncores)
        p0 = basic_p(bs$t0[1], bs$t[,1])
        #  old_bs_t <- c(old_bs_t, bs$t[,1])
      }
      #if(early_stop_R > 700){
      #    if(p0<0.01){
      # bs = boot::boot(df2,assoc, R = 2000, formula = formula,  stype = 'i', parallel = 'multicore', ncpus = ncores)
      # p0 = basic_p(bs$t0[1], bs$t[,1])
      #}
       #   }
    }
      #if(p0<0.001){
      # bs = boot::boot(df2,assoc, R = 5000, formula = formula, stype = ‘i’, parallel = “multicore”, ncpus = ncores)
      #  p0 = basic_p(bs$t0[1], bs$t[,1])
      #  }
      out <- data.frame(gene=gene,peak=this_peak,beta=coefs[1],se=coefs[2],z=coefs[3],p=coefs[4],boot_basic_p=p0)
      res<-rbind(res,out)
    }
  }
  #Update the SCENT.result field of the constructor in R:
  object@SCENT.result <- res
  return(object)
}
