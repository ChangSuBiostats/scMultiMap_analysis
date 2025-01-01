
# res_df_list is a list of two dataframe
# each of the dataframe should contain columns of peak, gene and p_var
get_overlap <- function(res_df_list, p_var='pval', p_cutoff = 0.05, subset_inds = NULL,
                       n_cutoff = NULL){
    grange_list <- sig_inds_list <- sig_grange_list <- list()
    if(!is.null(subset_inds)){
        for(i in 1:2){
            res_df_list[[i]] <- res_df_list[[i]][subset_inds,]
        }
    }
    for(i in 1:2){
        grange_list[[i]] <- Signac::StringToGRanges(res_df_list[[i]]$peak, sep = c('-', '-'))
        if(is.null(n_cutoff)){
            sig_inds_list[[i]] <- res_df_list[[i]][[p_var]] < p_cutoff
        }else{
            n_pair_inds <- 1:nrow(res_df_list[[i]])
            sig_inds_list[[i]] <- n_pair_inds %in% n_pair_inds[order(res_df_list[[i]][[p_var]], decreasing=F)[1:n_cutoff]]
        }
        sig_grange_list[[i]] <- grange_list[[i]][sig_inds_list[[i]]]
    }
    # -
    # evaluate overlap in peaks
    # -
    # background (all peaks) overlap
    bovp = with(GenomicRanges::findOverlaps(grange_list[[1]], grange_list[[2]]),
                   data.frame(d1_idx=queryHits,
                              d2_idx=subjectHits))
    # significant peak overlap
    sovp = with(GenomicRanges::findOverlaps(sig_grange_list[[1]], sig_grange_list[[2]]),
                   data.frame(d1_idx=queryHits,
                              d2_idx=subjectHits))
    # -
    # evaluate overlap in gene-peak pairs
    # -
    # background (all peak-gene pairs) overlap
    n_overlap_pairs <- sum(res_df_list[[1]]$gene[bovp$d1_idx] == 
                           res_df_list[[2]]$gene[bovp$d2_idx])
    # significant peak-gene pair overlap
    n_sig_both <- sum(res_df_list[[1]]$gene[sig_inds_list[[1]]][sovp$d1_idx] ==
                      res_df_list[[2]]$gene[sig_inds_list[[2]]][sovp$d2_idx])
    # overlap of (significant peak-gene pair in dataset 1) with (all pairs in dataset 2)
    n_sig_one_dataset <- numeric(2)
    for(i in 1:2){
        sovp_d1 = with(GenomicRanges::findOverlaps(sig_grange_list[[i]], grange_list[[2-i+1]]),
                       data.frame(d1_idx=queryHits,
                                  d2_idx=subjectHits))
        # significant & overlap with pairs in CG data
        n_sig_one_dataset[i] <- sum(res_df_list[[i]]$gene[sig_inds_list[[i]]][sovp_d1$d1_idx] == 
                                    res_df_list[[2-i+1]]$gene[sovp_d1$d2_idx])
    }
    # phyper: drawn white ball, white ball, black ball, drawn ball
    # drawn -- significant in dataset 2
    # white -- significant in dataset 1
    log_p_val <- phyper(n_sig_both, n_sig_one_dataset[1], 
                        n_overlap_pairs-n_sig_one_dataset[1], n_sig_one_dataset[2], 
                        lower.tail=F, log.p = T)
    # -n_sig_both
    # -n_sig_one_dataset[2]
    fc <- (n_sig_both/n_sig_one_dataset[1])/((n_sig_one_dataset[2])/(n_overlap_pairs))
    # (n_sig_both/n_sig_one_dataset[2])/((n_sig_one_dataset[1])/(n_overlap_pairs))
    return(list(log_p_val = log_p_val, 
                fc = fc,
                counts = list(n_sig_both, n_sig_one_dataset, n_overlap_pairs)))
}


load_proposed <- function(ct, datasets = c('PBMC_10k', 'PBMC_10k_nextgem', 'PBMC_10k_nextgem_X'),
                         fn_suffix = '', output_dir = './output'){
    if(length(ct) == 1) ct <- rep(ct, length(datasets))
    if(length(fn_suffix) == 1) fn_suffix <- rep(fn_suffix, length(datasets))
    names(ct) <- names(fn_suffix) <- datasets
    res_df_list <- list()
    for(dataset in datasets){
        res_df_list[[dataset]] <- read.table(sprintf("%s/%s/%s_proposed_observed_results%s.txt", output_dir, dataset, ct[dataset],
                                                    fn_suffix[dataset]))
        res_df_list[[dataset]]$pval <- 2 * pnorm(abs(res_df_list[[dataset]]$test_stat), lower.tail = F)
        res_df_list[[dataset]]$adj_pval <- p.adjust(res_df_list[[dataset]]$pval, method = 'BH')
    }
    return(res_df_list)
}

load_SCENT <- function(ct, n_sets = 10, datasets = c('PBMC_10k', 'PBMC_10k_nextgem', 'PBMC_10k_nextgem_X'),
                      fn_suffix = '', output_dir = './output'){
    if(length(ct) == 1) ct <- rep(ct, length(datasets))
    names(ct) <- datasets
    if(length(n_sets) == 1) n_sets <- rep(n_sets, length(datasets))
    names(n_sets) <- datasets
    if(length(fn_suffix) == 1) fn_suffix <- rep(fn_suffix, length(datasets))
    names(fn_suffix) <- datasets
    
    res_df_list <- list()
    for(dataset in datasets){
        if(fn_suffix[dataset] != '_permu_adjusted'){
            tmp <- list()
            for(i_set in 1:n_sets[dataset]){
                if(n_sets[dataset] != 50){
                    fn <- sprintf("%s/%s/%s_SCENT_n_sets_%i_i_set_%i_results%s.rds", 
                                            output_dir, dataset, ct[dataset], n_sets[dataset], i_set, fn_suffix[dataset])
                }else{
                    fn <- sprintf("%s/%s/%s_SCENT_i_set_%i_results%s.rds", 
                                                output_dir, dataset, ct[dataset], i_set, fn_suffix[dataset])
                }
                tmp[[i_set]] <- readRDS(fn)
            }
            res_df_list[[dataset]] <- do.call(rbind, tmp)
            res_df_list[[dataset]]$pval <- res_df_list[[dataset]]$boot_basic_p
            res_df_list[[dataset]]$adj_pval <- p.adjust(res_df_list[[dataset]]$pval, method = 'BH')
        }else{
            fn <- sprintf("%s/%s/%s_SCENT_observed_results%s.txt", 
                         output_dir, dataset, ct[dataset], fn_suffix[dataset])
            res_df_list[[dataset]] <- read.table(fn)
            res_df_list[[dataset]]$adj_zscore <- (res_df_list[[dataset]]$score - res_df_list[[dataset]]$bg.mean) / res_df_list[[dataset]]$bg.sd
            res_df_list[[dataset]]$pval <- 2 * pnorm(abs(res_df_list[[dataset]]$adj_zscore), lower.tail = F)
            res_df_list[[dataset]]$adj_pval <- p.adjust(res_df_list[[dataset]]$pval, method = 'BH')
        }
    }
    return(res_df_list)
}

load_Signac <- function(ct, datasets = c('PBMC_10k', 'PBMC_10k_nextgem', 'PBMC_10k_nextgem_X'), fn_suffix = '', output_dir = './output'){
    if(length(ct) == 1) ct <- rep(ct, length(datasets))
    names(ct) <- datasets
    if(length(fn_suffix) == 1) fn_suffix <- rep(fn_suffix, length(datasets))
    names(fn_suffix) <- datasets
    res_df_list <- list()
    for(dataset in datasets){
        res_df_list[[dataset]] <- read.table(sprintf("%s/%s/%s_Signac_observed_results%s.txt", output_dir, dataset, ct[dataset], fn_suffix[dataset]))
        res_df_list[[dataset]]$test_stat <- (res_df_list[[dataset]]$score - res_df_list[[dataset]]$bg.mean) / res_df_list[[dataset]]$bg.sd
        res_df_list[[dataset]]$pval <- 2 * pnorm(abs(res_df_list[[dataset]]$test_stat), lower.tail = F)
        res_df_list[[dataset]]$adj_pval <- p.adjust(res_df_list[[dataset]]$pval, method = 'BH')
    }
    return(res_df_list)
}

load_SEACells <- function(ct, datasets = c('PBMC_10k', 'PBMC_10k_nextgem', 'PBMC_10k_nextgem_X'), fn_suffix = '', output_dir = './output'){
    if(length(ct) == 1) ct <- rep(ct, length(datasets))
    names(ct) <- datasets
    res_df_list <- list()
    for(dataset in datasets){
        res_df_list[[dataset]] <- read.table(sprintf("%s/%s/%s_SEACells_observed_results.txt", output_dir, dataset, ct[dataset]), sep=',', header=T)
        print(sprintf('#NA SEACells p values: %i', sum(is.na(res_df_list[[dataset]]$pval))))
        res_df_list[[dataset]]$pval[is.na(res_df_list[[dataset]]$pval)] <- 1
        res_df_list[[dataset]]$adj_pval <- p.adjust(res_df_list[[dataset]]$pval, method = 'BH')
    }
    return(res_df_list)
}

binarize_importance <- function(gbm_res, binarize_method = 'top_n', method_par = 5){
    genes <- unique(gbm_res$gene)
    pos_inds <- numeric()
    for(g in genes){
        g_index <- gbm_res$g_index[gbm_res$gene == g]
        if(binarize_method == 'top_n'){
            pos_inds <- c(pos_inds, g_index[1:min(method_par, length(g_index))])
        }else if(binarize_method == 'quantile'){
            pos_inds <- c(pos_inds, g_index[1:((1-method_par)*length(g_index))])
        }
    }
    gbm_res$pval <- 1
    gbm_res$pval[pos_inds] <- 0
    return(gbm_res)
}

load_GBM <- function(ct, datasets, binarize_method = 'quantile', method_par = 0.95, output_dir = './output'){
    res_df_list <- list()
    for(dataset in datasets){
        res_df_list[[dataset]] <- readRDS(sprintf('%s/%s/GBM_%s.rds', output_dir, dataset, ct))
        res_df_list[[dataset]]$g_index <- 1:nrow(res_df_list[[dataset]])
        res_df_list[[dataset]] <- binarize_importance(res_df_list[[dataset]], binarize_method, method_par)
    }
    return(res_df_list)
}

single_ct_plot <- function(df, title){
    g <- ggplot(df) + 
    geom_bar(aes(x = setting, y = x, fill = method), stat = 'identity', position = 'dodge') +
    labs(title = title, y = title) + 
    theme_classic(base_size = 20)
    return(g)
}

multi_ct_plot <- function(df, title){
    g <- ggplot(df) + 
    geom_bar(aes(x = ct, y = x, fill = method), stat = 'identity', position = 'dodge') +
    labs(title = title, y = title) + 
    theme_classic(base_size = 20)
    if(title == 'fc') g <- g +  coord_cartesian(ylim=c(0.98,2.2))
    return(g)
}