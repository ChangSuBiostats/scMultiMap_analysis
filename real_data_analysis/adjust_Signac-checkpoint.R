# adjust for batch effects in Signac estimates
# by constructing a null distribution that captures the confounding effects of batches on Signac estimates;
#  ``` Rscript run_Signac.R --permu_within_batch=TRUE --i_permu=1 ... ```
# use the null distribution to derive a new p-value for the observed Signac estimates via z-score transformation;
# and use it in downstream analyses such as reproducibility with plac-seq data
library(dplyr)
library("optparse")
option_list = list(
    make_option(c("--study"), type="character", default="brain_CG",
              help="which dataset to analyze", metavar="character"),
    make_option(c("--i_ct"), type="integer", default=1,
              help="which cell type", metavar="integer"),
    make_option(c("--cell_subset"), type="character", default="Control",
                help="which subset of cells to analyze", metavar="character"),
    make_option(c("--method"), type="character", default="Signac",
              help="which method to analyze", metavar="character"),
    make_option(c("--run_TF_peak"), type="character", default="FALSE",
                help="whether to run TF peak association")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
study <- opt$study
i_ct <- opt$i_ct
cell_subset <- opt$cell_subset
method <- opt$method
fn_suffix <- ifelse(cell_subset == 'Control', '', '_All_cells')
TF_peak_suffix <- ifelse(opt$run_TF_peak == 'TRUE', '_TF_peak', '')
print(cell_subset)
print(fn_suffix)
print(method)
print(TF_peak_suffix)

# load data
if(grepl('brain', study)){
    pn <- 'peaks'
    if(study == 'brain_CG'){
        cts <- c('Microglia', 'Excitatory', 'Inhibitory', 'Oligodendrocytes', 'Astrocytes')
    }else if(study == 'brain_ROSMAP'){
        cts <- c('Mic', 'Exc', 'Inh', 'Oli', 'Ast')
    }
    if(cell_subset == 'Control'){
        obj_suffix <- '_control'
    }else if(cell_subset == 'All'){
        obj_suffix <- ''
    }
}else if(grepl('PBMC', study)){
    pn <- 'peaks_ct'
    obj_suffix <- ''
    cts <- c('CD14 Mono', 'CD4 TCM', 'CD8 Naive', 'CD4 Naive', 'CD8 TEM')
}
ct <- cts[i_ct]
print(ct)
output_dir <- sprintf('output/%s/%s', study, ct)

if(method == 'Signac'){
obs_res <- sprintf('%s_Signac_observed_results%s%s.txt', output_dir, fn_suffix, TF_peak_suffix) %>% read.table
permu_res <- array(dim = c(100, nrow(obs_res), 3))
for(i_permu in 1:100){
    tmp <- sprintf('%s_Signac_observed_results_permu_within_batch_i_%s%s%s.txt', output_dir, i_permu, fn_suffix, TF_peak_suffix) %>% read.table 
    permu_res[i_permu,,] <- tmp[,1:3] %>% as.matrix
}
# Signac estimates on observed data
obs_zscore <- (obs_res[[3]] - obs_res[[1]])/ obs_res[[2]]
# Signac estimates on on permuted data
# where peak counts were randomly permuted within each batch (sample)
# such that it provides a null distribution with where peak & gene are independent, while preserving the confounding effects of batch on Signac zscores
permu_zscore <- t((permu_res[,,3]-permu_res[,,1])/permu_res[,,2])

}else if(method == 'SCENT'){
  print('load SCENT results')
  source('reproducibility_helper.R')
  obs_res <- load_SCENT('Microglia', 50, 'brain_CG')[[1]]
  permu_zscore <- matrix(nrow = 100, ncol = nrow(obs_res))
  for(i_permu in 1:100){
    tmp <- load_SCENT('Microglia', 10, 'brain_CG', sprintf('_permu_within_batch_i_%i', i_permu))[[1]]
    permu_zscore[i_permu,] <- tmp$z
      permu_zscore_bootstrap[i_permu,] <- qnorm(tmp$boot_basic_p)
  }
  obs_zscore <- obs_res$z
  permu_zscore <- t(permu_zscore)
    permu_zscore_bootstrap <- t(permu_zscore_bootstrap)
}

tmp_df <- data.frame(bg.mean = rowMeans(permu_zscore), 
                     bg.sd = apply(permu_zscore, 1, sd),
                    score = obs_zscore,
                    peak = obs_res$peak,
                    gene = obs_res$gene)
if(method == 'SCENT'){
    tmp_df$boot_basic_p = obs_res$boot_basic_p
    permu_zscore_bootstrap[is.infinite(permu_zscore_bootstrap)] <- -5 # avoid pvalue=0, zscore=-Inf
    tmp_df$permu_zscore_bootstrap.bg.mean = rowMeans(permu_zscore_bootstrap)
    tmp_df$permu_zscore_bootstrap.bg.sd = apply(permu_zscore_bootstrap, 1, sd)
}

adj_fn <- sprintf('%s_%s_observed_results_permu_adjusted%s%s.txt', output_dir, method, fn_suffix, TF_peak_suffix)
write.table(tmp_df, adj_fn, quote = F)
print(adj_fn)
