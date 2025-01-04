TF_motif_enrich <- function(peaks, sparsity, obj){
    enriched.motifs <- FindMotifs(
      object = obj,
      features = peaks
    )
    enriched.motifs$sparsity <- sparsity[match(enriched.motifs$motif.name, names(sparsity))]
    return(enriched.motifs)
}

make_trio <- function(TF, results, TF_target_p, TF_target_coexp,
                      TF_peak_p_df, 
                      pair_sigp_cutoff = 0.05, simul = T,
                     subset = T, eval_proposed=T){
    #results$coexp_gene_peak <- results$covar / sqrt(results$sigma_sq_gene * results$sigma_sq_peak)
    # extract all peaks that contain the motif for this TF
    peaks_w_binding_motifs <- names(which(motif.object@data[,full_TF_names == TF] != 0))
    results_peaks_w_motif_inds <- results$peak %in% peaks_w_binding_motifs
    if(any(results_peaks_w_motif_inds)){
        # extract three p values, focusing on peaks and genes related to the TF
        # from peak to gene
        peak_gene_pval <- results$pval[results_peaks_w_motif_inds]
        # from gene to TF
        gene_TF_pval <- TF_target_p[TF, results$gene][results_peaks_w_motif_inds]
        # calculate the associations between TF and peak
        peak_TF_pval <- TF_peak_p_df$pval[match(results$peak[results_peaks_w_motif_inds], rownames(TF_peak_p_df))]

        if(eval_proposed){
            # extract co-expressions
        # from peak to gene
        peak_gene_coexp <- results$coexp[results_peaks_w_motif_inds]
        # from gene to TF
        gene_TF_coexp <- TF_target_coexp[TF, results$gene][results_peaks_w_motif_inds]
        # from peak to TF
        peak_TF_coexp <- TF_peak_p_df$coexp[match(results$peak[results_peaks_w_motif_inds], rownames(TF_peak_p_df))]
        }
        
        # criterion for a trio
        if(simul){
            trio_inds <- (peak_gene_pval < pair_sigp_cutoff) & (gene_TF_pval < pair_sigp_cutoff) & (peak_TF_pval < pair_sigp_cutoff)
        }else{
            trio_inds <- sapply((peak_gene_pval < pair_sigp_cutoff) + (gene_TF_pval < pair_sigp_cutoff) + (peak_TF_pval < pair_sigp_cutoff),
                               function(x) x >= 2)
        }
        if(any(trio_inds)){
            tmp <- results[results_peaks_w_motif_inds, c('peak','gene','pval')]
            colnames(tmp)[colnames(tmp) == 'pval'] <- 'pval_peak_gene'
            tmp$pval_TF_gene <- gene_TF_pval
            tmp$pval_peak_TF <- peak_TF_pval
            if(eval_proposed){
            tmp$coexp_peak_gene <- peak_gene_coexp
            tmp$coexp_TF_gene <- gene_TF_coexp
            tmp$coexp_peak_TF <- peak_TF_coexp
            }
            tmp$TF <- TF
        }
        if(!any(trio_inds)){
            return()
        }else if(subset){
            return(tmp[trio_inds,])
        }else{
            return(tmp)
        }
    }
}


gene_GO <- function(genes, bg=NULL, source=NULL, as_short_link = FALSE){
    # https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis
    gostres <- gost(query = genes,
                #proposed$gene[abs(proposed$test_stat) > 1.96] %>% unique, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = bg, #top_genes[[1]], #NULL, 
                numeric_ns = "", sources =source, as_short_link = as_short_link, highlight = TRUE)
    if(!as_short_link){
        return(gostres$result)
    }else{
        print(gostres)
    }
}

reorder_trios_by_TF <- function(coexp_mat, top_TFs, full_TFs,
                               coexp_cutoff = 0.1){
    sub_inds <- (abs(coexp_mat[,'coexp_TF_gene']) > coexp_cutoff) & 
                                (abs(coexp_mat[,'coexp_peak_gene']) > coexp_cutoff) &
                                (abs(coexp_mat[,'coexp_peak_TF']) > coexp_cutoff)
    full_TFs <- full_TFs[sub_inds]
    sub_coexp_mat <- coexp_mat[sub_inds,
                              c('coexp_peak_gene',
                                  'coexp_peak_TF',
                                  'coexp_TF_gene')]
    for(tf in top_TFs){
        otter_dendro <- hclust(d = dist(x = sub_coexp_mat[full_TFs==tf,]))
        sub_coexp_mat[full_TFs==tf, ] <- sub_coexp_mat[full_TFs==tf, ][otter_dendro$order,]
    }
    return(list(sub_coexp_mat, full_TFs))
}

reorder_trios <- function(coexp_mat, top_TFs, full_TFs,
                               coexp_cutoff = 0.1){
    sub_inds <- (abs(coexp_mat[,'coexp_TF_gene']) > coexp_cutoff) & 
                                (abs(coexp_mat[,'coexp_peak_gene']) > coexp_cutoff) &
                                (abs(coexp_mat[,'coexp_peak_TF']) > coexp_cutoff)
    full_TFs <- full_TFs[sub_inds]
    sub_coexp_mat <- coexp_mat[sub_inds,
                              c('coexp_peak_gene',
                                  'coexp_peak_TF',
                                  'coexp_TF_gene')]
    #for(tf in top_TFs){
        otter_dendro <- hclust(d = dist(x = sub_coexp_mat))
        sub_coexp_mat <- sub_coexp_mat[otter_dendro$order,]
    #}
    return(list(sub_coexp_mat, full_TFs[otter_dendro$order]))
}

plot_heatmap <- function(coexp_mat, full_TFs, top_TFs, ct){
    # Annotations ===================================================
    # create a data frame for column annotation
    ann_df <- data.frame(Group = c("peak-gene", "peak-TF", "TF-gene"))
    row.names(ann_df) <- colnames(coexp_mat)

    gene_functions_df <- data.frame(TF = full_TFs)
    row.names(gene_functions_df) <- rownames(coexp_mat)
    ann_colors <- list(
        # https://coolors.co/palettes/popular/6%20colors
      gene_functions = c("#ef476f",
                         "#f78c6b",
                         "#ffd166",
                         "#06d6a0", 
                         "#118ab2",
                         "#073b4c"), 
        Group = c("peak-gene" = "#26547c",
                  "peak-TF" = "#ef476f",
                  "TF-gene" = "#ffd166"))
    names(ann_colors[[1]]) <- top_TFs

    pheatmap(coexp_mat,
                      col = rev(brewer.pal(11, 'RdBu')), # choose a colour scale for your data
                      cluster_rows = F, cluster_cols = F, # set to FALSE if you want to remove the dendograms
                      #clustering_distance_cols = 'euclidean',
                      #clustering_distance_rows = 'euclidean',
                      #clustering_method = 'ward.D',
                      annotation_row = gene_functions_df, # row (gene) annotations
                      annotation_col = ann_df, # column (sample) annotations
                      annotation_colors = ann_colors, # colours for your annotations
                      annotation_names_row = F, 
                      annotation_names_col = F,
                      fontsize_row = 10,          # row label font size
                      fontsize_col = 7,          # column label font size 
                      angle_col = 45, # sample names at an angle
                      legend_breaks = seq(-1,1,0.2), #c(-1, 0, 1), # legend customisation
                      #legend_labels = c("Low", "Medium", "High"), # legend customisation
                      show_colnames = T, show_rownames = F, # displaying column and row names
                      main = ct) # a title for our heatmap

}
                                
gene_GO <- function(genes, bg=NULL, source=NULL, as_short_link = FALSE){
    # https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis
    gostres <- gost(query = genes,
                #proposed$gene[abs(proposed$test_stat) > 1.96] %>% unique, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = bg, #top_genes[[1]], #NULL, 
                numeric_ns = "", sources =source, as_short_link = as_short_link, highlight = TRUE)
    if(!as_short_link){
        return(gostres$result)
    }else{
        print(gostres)
    }
}

load_trios_list <- function(p_top_gene, ct, study, trio_pval_cutoff, output_dir){
    pairs_suffix <- ifelse(p_top_gene == 5000, '_p_top_5000_50000', '')
    fn <- sprintf('%s/%s/trios_%s%s%s.rds', 
                  output_dir,
                  study,
                  ct, 
                   pairs_suffix,
                   sprintf('_pval_cutoff_%.2e', trio_pval_cutoff))
    #print(fn)
    TF_trios_list <- readRDS(fn)
    trios_list <- TF_trios_list[[1]]
    return(trios_list)
}
                                
summarize_enr <- function(trios_list, p_top = 5, p_TF_top = 10, verbose = F){
    gene_enrich_short_list <- gene_enrich_list <- list()
    for(TF in names(table(trios_list$TF) %>% sort(decreasing = T))[1:min(p_TF_top, length(unique(trios_list$TF)))]){
        print(TF)
        gene_enrich_list[[TF]] <- gene_GO(trios_list$gene[trios_list$TF == TF] %>% unique,
                                          source = c('GO:BP', 'GO:MF', 'GO:CC'))
        if(verbose) gene_GO(trios_list$gene[trios_list$TF == TF] %>% unique,
                                          source = c('GO:BP', 'GO:MF', 'GO:CC'), as_short_link=T)
        # save only the 'highlighted' / driver GO terms
        gene_enrich_list[[TF]] <- gene_enrich_list[[TF]][gene_enrich_list[[TF]]$highlighted & gene_enrich_list[[TF]]$source == 'GO:BP', ]
        if(length( gene_enrich_list[[TF]]) > 0){
            # focus on top ones
            tmp <- gene_enrich_list[[TF]][order(gene_enrich_list[[TF]]$p_value, decreasing = F),]
            gene_enrich_short_list[[TF]] <- tmp[1:min(p_top, nrow(tmp)),]
            gene_enrich_short_list[[TF]]$TF <- TF
        }
        
    }
    gene_enrich_short_df <- do.call(rbind, gene_enrich_short_list)
    gene_enrich_short_df <- gene_enrich_short_df[!is.na(gene_enrich_short_df[,1]),]
    if(verbose) table(gene_enrich_short_df$term_name) %>% sort(decreasing = T) %>% print
    return(gene_enrich_short_df)
}
                                    
make_pheatmap <- function(df){
    all_terms <- unique(df$term_name)
    enrich_mat <- matrix(0, nrow = length(all_terms), ncol = length(unique(df$TF)))
    rownames(enrich_mat) <- all_terms
    colnames(enrich_mat) <- unique(df$TF)
    for(i in 1:nrow(df)){
        enrich_mat[df$term_name[i] == all_terms, 
               colnames(enrich_mat) == df$TF[i]] <- -log10(df$p_value[i])
    }
    return(enrich_mat)
}