load.plac <- function(igr1, igr2, egr, peaks_ranges, promoters_ranges, verbose=F) {
    # reference: https://github.com/KellisLab/AD_regulome_analysis/blob/main/integration_linking_modules/fit_adgwas.R
    #idf = read.table(interactome_file, sep="\t", header=TRUE)
    #igr1 = with(idf, GRanges(seqnames=chr1, ranges=IRanges(start1, end1)))
    #igr2 = with(idf, GRanges(seqnames=chr2, ranges=IRanges(start2, end2)))

    #egr = with(read.table(enhancers_file, sep="\t"), GRanges(seqnames=V1, ranges=IRanges(V2, V3)))
    eovp = findOverlaps(egr, peaks_ranges)
    egr = peaks_ranges[unique(subjectHits(eovp))] ### subset peaks to overlapping enhancers (comment by the author)
    if(verbose){
        print(length(peaks_ranges))
        print(length(egr))
    }
    # process interactome file by only keeping the promotor-enhancer link
    eovp_12 = with(findOverlaps(egr, igr1),
                   data.frame(enh_idx=queryHits,
                              idf_idx=subjectHits))
    eovp_21 = with(findOverlaps(egr, igr2),
                   data.frame(enh_idx=queryHits,
                              idf_idx=subjectHits))
    govp_21 = with(findOverlaps(promoters_ranges, igr1),
                   data.frame(prom_idx=queryHits,
                              idf_idx=subjectHits))
    govp_12 = with(findOverlaps(promoters_ranges, igr2),
                   data.frame(prom_idx=queryHits,
                              idf_idx=subjectHits))
    eg = rbind(merge(eovp_12, govp_12)[c("prom_idx", "idf_idx", "enh_idx")],
               merge(eovp_21, govp_21)[c("prom_idx", "idf_idx", "enh_idx")])
    if(verbose){
        print(dim(eg))
        print(length(igr1))
    }
    d2n = distanceToNearest(egr[eg$enh_idx], promoters_ranges[eg$prom_idx])
    eg_distance = distance(egr[eg$enh_idx], promoters_ranges[eg$prom_idx])
    closest_gene_distance = as.data.frame(d2n)$distance
    plac = with(eg, data.frame(gene=names(promoters_ranges)[prom_idx],
                               fdr=igr1$fdr[idf_idx],
                               ClusterSummit=igr1$ClusterSummit[idf_idx],
                               ClusterNegLog10P=igr1$ClusterNegLog10P[idf_idx],
                               ClusterLabel=igr1$ClusterLabel[idf_idx],
                               ClusterSize=igr1$ClusterSize[idf_idx],
                               distance=eg_distance,
                               is_closest_gene=closest_gene_distance >= eg_distance,
                               peak=names(egr)[enh_idx]))
    plac = plac[order(plac$fdr),]
    plac = plac[!duplicated(plac[c("peak", "gene")]),]
    rownames(plac) = paste0(plac$peak, "|", plac$gene)
    return(plac)
}

grange_to_str <- function(gr){
    paste(seqnames(gr), ranges(gr)%>% start, ranges(gr)%>% end, sep = '-')
}

make_grange_hg19_to_38 <- function(df, seqf = 'chr', startf = 'start', endf = 'end', to_reduce=F, verbose=F){
    grange_hg19 <- GenomicRanges::makeGRangesFromDataFrame(df = df, 
                                                         seqnames.field = seqf,
                                                         start.field = startf,
                                                         end.field = endf,
                                                          keep.extra.columns=T)

    grange_hg19$id <- 1:length(grange_hg19)
    seqlevelsStyle(grange_hg19) = "UCSC"  # necessary
    grange_hg38 <- liftOver(grange_hg19, chain)
    if(to_reduce){
        for(i in 1:length(grange_hg38)){
            # reduce first orders the ranges in x from left to right, then merges the overlapping or adjacent ones.
            # https://www.biostars.org/p/464758/
            if(length(grange_hg38[[i]]) > 1) grange_hg38[[i]] <- reduce(grange_hg38[[i]])
        }
    }else{
        grange_hg38 <- unlist(grange_hg38)
    }
    genome(grange_hg38) = "hg38"#"GRCh38"
    if(verbose){
        length(grange_hg38) %>% print
        length(grange_hg19) %>% print
        print(sprintf('ratio of hg38 matched peaks versus hg19 peaks: %.2f', length(grange_hg38)/length(grange_hg19)))
    }
    return(grange_hg38)
}

# validation with HiChIP or PCHiC data
validate_with_assay <- function(res_list, i_dataset, cts, assay = 'HiChIP', pvar = 'adj_pval', p_cutoff = 0.2,
                                assay_sig_cutoff = list(HiChIP = 10^(-3), PCHiC = 5),
                               dist_cutoff = NULL, text_comp = F){
    pval <- n_overlap <- enr <- ratio <- oddsratio <-  numeric()
    overlap_pairs <- list()
    for(ct in cts){
        #print(ct)
        overlap_pairs_ct <- list()
        for(method in names(res_list[[ct]])){
            res_df <- res_list[[ct]][[method]][[i_dataset]]
            if(is.null(dist_cutoff)){
                dist_subset <- rep(T, nrow(res_df))
            }else{
                dist_subset <- selected_pair_list[[ct]][[i_dataset]]$gene_peak_dist <= dist_cutoff
                print(sprintf('subset to %i peak-gene pairs within %e base pairs', sum(dist_subset), dist_subset))
            }
            res_df <- res_df[dist_subset,]
            grange_peaks <- Signac::StringToGRanges(res_list[[ct]][[1]][[i_dataset]]$peak[dist_subset], sep = c('-', '-'))
            sig_res_inds <- res_df[[pvar]] < p_cutoff
            if(assay == 'HiChIP'){
                sig_assay_df <- interacting_peaks_hg38[interacting_peaks_hg38@elementMetadata[[paste0('qval_', cts_str[ct])]] < assay_sig_cutoff[[assay]]]
            }else if(assay == 'PCHiC'){
                # filter for interactions with CHiCAGO scores >=5 in the cell type of interests
                sig_assay_df <- PCHiC_peaks_hg38[PCHiC_peaks_hg38@elementMetadata[[PCHiC_cts[ct]]] > assay_sig_cutoff[[assay]]]
            }else if(assay == 'eQTL'){
                # consider only the 'bg_genes', i.e. genes considered in results / peak-gene association analysis
                # and focus on cis-SNPs significant in Naive CD14 Monocytes
                # which is more consistent with the biology of 10X PBMC sample from healthy subjects
                eQTL_tab <- eQTL_tab[(eQTL_tab$Gene %in% res_df$gene) & !is.na(eQTL_tab$Naive.p.value),]
                new_eQTL <- data.frame(chr = paste0('chr', eQTL_tab$SNP.Chrm), snp.pos = eQTL_tab$SNP.pos, gene = eQTL_tab$Gene)
                sig_assay_df <- make_grange_hg19_to_38(new_eQTL, "chr", "snp.pos", "snp.pos")
            }
                
            bovp = with(GenomicRanges::findOverlaps(grange_peaks, sig_assay_df),
                           data.frame(d1_idx=queryHits,
                                      d2_idx=subjectHits))
            sovp = with(GenomicRanges::findOverlaps(grange_peaks[sig_res_inds], 
                                                sig_assay_df),
                           data.frame(d1_idx=queryHits,
                                      d2_idx=subjectHits))

            if(assay == 'HiChIP'){
                # number of HiChIP pairs among significant peak-gene pairs
                assay_sig_inds <- res_df$gene[sig_res_inds][sovp$d1_idx] == sig_assay_df$GeneName[sovp$d2_idx]
                n_assay_sig <- sum(assay_sig_inds)
                # number of HiChIP pairs among all tested peak-gene pairs
                n_assay <- sum(res_df$gene[bovp$d1_idx] == sig_assay_df$GeneName[bovp$d2_idx])
                # consider peak-gene pairs where genes overlapped with the HiChIP data
                res_gene_inds <- (res_df$gene %in% sig_assay_df$GeneName)
            }else if(assay == 'PCHiC'){
                if(text_comp){
                    if(sum(sovp$d2_idx) != 0){
                        sbait <- strsplit(sig_assay_df$baitName[sovp$d2_idx], ';')
                        assay_sig_inds <- sapply(which(!is.na(sbait)), function(i) any(res_df$gene[sig_res_inds][sovp$d1_idx][i] == sbait[[i]]))
                                                 print(assay_sig_inds)
                        n_assay_sig <- sum(assay_sig_inds, na.rm=T)
                    }else{
                        n_assay_sig <- 0
                    }
                    sbait <- strsplit(sig_assay_df$baitName[bovp$d2_idx], ';')
                    n_assay <- sapply(which(!is.na(sbait)), function(i) any(res_df$gene[bovp$d1_idx][i] == sbait[[i]])) %>% sum(na.rm=T)
                                      print(n_assay)
                }else{
                    assay_sig_inds <- res_df$gene[sig_res_inds][sovp$d1_idx] == sig_assay_df$baitName[sovp$d2_idx]
                    n_assay_sig <- sum(assay_sig_inds, na.rm=T)
                    n_assay <- sum(res_df$gene[bovp$d1_idx] == sig_assay_df$baitName[bovp$d2_idx], na.rm=T)
                }
                res_gene_inds <- (res_df$gene %in% sig_assay_df$baitName)
                #res_gene_inds <- rep(T, length(res_df$gene))
            }else if(assay == 'eQTL'){
                assay_sig_inds <- res_df$gene[sig_res_inds][sovp$d1_idx] == sig_assay_df$gene[sovp$d2_idx]
                n_assay_sig <- sum(assay_sig_inds, na.rm=T)
                n_assay <- sum(res_df$gene[bovp$d1_idx] == sig_assay_df$gene[bovp$d2_idx], na.rm=T)
                res_gene_inds <- (res_df$gene %in% sig_assay_df$gene)
            }
            n_not_assay <- nrow(res_df[res_gene_inds,]) - n_assay
            n_sig <- sum(sig_res_inds & res_gene_inds)
            #n_not_assay <- nrow(res_df) - n_assay # length(unique(bovp$d1_idx)) - n_assay
            #n_sig <- sum(sig_res_inds) # sum(sig_res_inds[unique(bovp$d1_idx)])

            #print(c(n_assay_sig, n_assay, n_not_assay, n_sig))
            #phyper(n_assay_sig, n_assay, n_not_assay, n_sig, lower.tail=F) %>% print
        
            pval <- c(pval, phyper(n_assay_sig, n_assay, n_not_assay, n_sig, lower.tail=F))
            n_overlap <- c(n_overlap, n_assay_sig)
            ratio <- c(ratio, (n_assay_sig/n_sig))
            enr <- c(enr, (n_assay_sig/n_sig) / (n_assay/(n_assay+n_not_assay)))
                                     # print(n_assay/(n_assay+n_not_assay))
            A <- n_assay_sig
            B <- n_sig - A
            C <- n_assay - A
            D <- (n_assay+n_not_assay) - (A + B + C)
            oddsratio <- c(oddsratio, (A*D)/(B*C))
            if(length(sovp$d2_idx) == 0){
                overlap_pairs_ct[[method]] <- paste(res_df$gene, res_df$peak, sep='|')[sig_res_inds][sovp$d1_idx][assay_sig_inds]
            }else{
                overlap_pairs_ct[[method]] <-  ''
            }
        }
        pval[n_overlap == 0] <- 1
        overlap_pairs[[ct]] <- overlap_pairs_ct
    }
    return(list(pval = pval, n_overlap = n_overlap, ratio = ratio, enr = enr, n_assay = n_assay,
                oddsratio = oddsratio,
               overlap_pairs = overlap_pairs))
}

make_plot <- function(x, title, cts, methods){
    x_df <- data.frame(x = x, 
                       ct = rep(cts, each = length(methods)),
                       method = rep(methods, length(cts)))
    x_df$ct <- factor(x_df$ct, levels = cts, labels = cts)
    g <- ggplot(x_df) +
        geom_bar(aes(y = x, x = method), stat='identity', position = 'dodge', show.legend = FALSE) +
        facet_wrap(~ct, ncol = length(cts), scales = "free") +
        theme_classic(base_size=20) +
        labs(title = title)
    print(g)
    return(g)
}

validate_with_assay_plots <- function(res, cts, methods){
    require("ggVennDiagram")
    g_list <- list()
    for(ct in cts){
        g_list[[ct]] <- ggVennDiagram(res$overlap_pairs[[ct]]) + labs(title = ct) +
        scale_fill_gradient(low = "white", high = "red")    
    }
    options(repr.plot.width = 4*length(cts), repr.plot.height = 4)
    grid.arrange(grobs = g_list, nrow = 1)

    make_plot(-log10(res$pval), '-log p-value', cts, methods)
    make_plot(res$n_overlap, 'number of identified pairs',  cts, methods)
    make_plot(res$enr, 'Enrichment',  cts, methods)
}

count_pair_match <- function(gr1, gr2, gene1, gene2){
    matching <- findOverlaps(gr1, gr2)
    gene_match_bin <- gene1[queryHits(matching)] == gene2[subjectHits(matching)]
    # number of unique peaks that have their genes matched with ground truth in gr1
    return(list(c(n_matched_peaks = length(queryHits(matching)),
            n_matched_peak_genes_pairs = length((queryHits(matching))[gene_match_bin])),
               queryHits(matching)[gene_match_bin],
               gr2[subjectHits(matching)][gene_match_bin]))
}
                                      
validate_with_SNP_gene_pairs <- function(dataset, results, pvar = 'adj_pval', p_cutoff = 0.2, verbose = F,
                                        split_gene_set = F, take_max = F){
    require(openxlsx)
    if(dataset == 'Science.CD14.Mono'){
        # http://dx.doi.org/10.1126/science.1246949
        # B. P. Fairfax et al., Science 343, 1246949 (2014). DOI: 10.1126/science.1246949
        # eQTL in CD14 Monocytes
        eQTL_tab <- read.xlsx('data/validation//PBMC//1246949stables2.xlsx', sheet = 'B. Complete.cis.summary.all.ind') # A. cis.peak.eSNPs.228.paired
        # consider only the 'bg_genes', i.e. genes considered in results / peak-gene association analysis
        # and focus on cis-SNPs significant in Naive CD14 Monocytes
        # which is more consistent with the biology of 10X PBMC sample from healthy subjects
        eQTL_tab <- eQTL_tab[(eQTL_tab$Gene %in% results$gene) & !is.na(eQTL_tab$Naive.p.value),]
        if(verbose){
            print(summary(eQTL_tab$Naive.p.value)) 
            print(summary(eQTL_tab$Naive.FDR)) # all gene-SNP pairs in this table are significant in Naive CD14 Monocytes, i.e. representing cis-eQTLs 
        }
        SNP_gene_pair <- data.frame(chr = paste0('chr', eQTL_tab$SNP.Chrm), snp.pos = eQTL_tab$SNP.pos, gene = eQTL_tab$Gene)
    }else if(dataset == 'NG.2023.Microglia.GWAS.pcHiC'){
        # priortized GWAS variants and their target genes validated by pcHiC
        pg_tab <- read.xlsx('/projects/chang/scGRN/peak_gene/data/brain/microglia_NG_2023/41588_2023_1506_MOESM3_ESM.xlsx',
                           sheet = 'Supplementary Table 6d', startRow = 2)
        #pg_tab <- pg_tab[pg_tab[['Gene(s).whose.promoter.falls.in.the.interacting.pcHi-C.bin']] %in% results$gene,]
        
        if(split_gene_set){
            pg_tab_list <- list()
            for(i in 1:nrow(pg_tab)){
                tmp_genes <- strsplit(pg_tab[['Gene(s).whose.promoter.falls.in.the.interacting.pcHi-C.bin']][i], ';')[[1]]
                
                tmp_df <- data.frame(gene = tmp_genes,
                                     chr = paste0('chr', pg_tab$chr[i]),
                                     snp.pos = pg_tab[['pos.hg38']][i],
                                    rsid = pg_tab[['rsid']][i])
                pg_tab_list[[i]] <- tmp_df
            }
            SNP_gene_pair <- do.call(rbind, pg_tab_list)
        }else{
            SNP_gene_pair <- data.frame(chr = paste0('chr', pg_tab$chr), 
                                 snp.pos = pg_tab[['pos.hg38']],
                                gene = pg_tab[['Gene(s).whose.promoter.falls.in.the.interacting.pcHi-C.bin']])
        }
        # there could be multiple Hi-C links linked to the same gene
        SNP_gene_pair <- unique(SNP_gene_pair)
    }else if(dataset == 'NG.2023.Microglia.GWAS.pcHiC.peak'){
        # the peaks that overlap with priortized GWAS variants and their target genes validated by pcHiC
        pg_tab <- read.xlsx('/projects/chang/scGRN/peak_gene/data/brain/microglia_NG_2023/41588_2023_1506_MOESM3_ESM.xlsx',
                           sheet = 'Supplementary Table 6d', startRow = 2)
        # expand sets of genes into seperate gene-peak pairs
        
        if(split_gene_set){
            pg_tab_list <- list()
            for(i in 1:nrow(pg_tab)){
                tmp_genes <- strsplit(pg_tab[['Gene(s).whose.promoter.falls.in.the.interacting.pcHi-C.bin']][i], ';')[[1]]
                
                tmp_df <- data.frame(gene = tmp_genes,
                                     peak.start = pg_tab[['Residing.pcHi-C.bin.start.pos.(hg19)']][i],
                                     peak.end = pg_tab[['Residing.pcHi-C.bin.end.pos.(hg19)']][i],
                                     chr = paste0('chr', pg_tab$chr[i]),
                                     snp = pg_tab$rsid[i],
                                    gwas_locus = pg_tab[['GWAS.locus']][i])
                pg_tab_list[[i]] <- tmp_df
            }
            SNP_gene_pair <- do.call(rbind, pg_tab_list)
        }else{
            #pg_tab <- pg_tab[pg_tab[['Gene(s).whose.promoter.falls.in.the.interacting.pcHi-C.bin']] %in% results$gene,]
            SNP_gene_pair <- data.frame(chr = paste0('chr', pg_tab$chr), 
                                    peak.start = pg_tab[['Residing.pcHi-C.bin.start.pos.(hg19)']],
                                    peak.end = pg_tab[['Residing.pcHi-C.bin.end.pos.(hg19)']], 
                                    gene = pg_tab[['Gene(s).whose.promoter.falls.in.the.interacting.pcHi-C.bin']],
                                    snp = pg_tab$rsid,
                                       gwas_locus = pg_tab[['GWAS.locus']])
        }
    }else if(dataset == 'NG.2023.Microglia.GWAS.priortized'){
        # priortized GWAS variants and their target genes validated by pcHiC
        pg_tab <- read.xlsx('/projects/chang/scGRN/peak_gene/data/brain/microglia_NG_2023/41588_2023_1506_MOESM3_ESM.xlsx',
                           sheet = 'Supplementary Table 6c', startRow = 2)
        # expand sets of genes into seperate gene-peak pairs
        eqtl_cols <- which(grepl('eQTL', colnames(pg_tab)))[1:3]
        pg_tab_list <- list()
        for(i in 1:nrow(pg_tab)){
            tmp_genes <- lapply(eqtl_cols, function(j) strsplit(pg_tab[[j]][i],';')[[1]]) %>% unlist %>% unique
            tmp_genes <- tmp_genes[!is.na(tmp_genes)]
            if(length(tmp_genes)==0) next()
            if(take_max){
                tmp_genes <- lapply(eqtl_cols, function(j) strsplit(pg_tab[[j]][i],';')[[1]]) %>% unlist 
                tmp_genes <- names(table(tmp_genes))[which.max(table(tmp_genes))]
            }             
            tmp_df <- data.frame(gene = tmp_genes,
                                snp.pos = pg_tab[['pos.hg38']][i],
                                chr = paste0('chr', pg_tab$chr[i]),
                                    snp = pg_tab$rsid[i],
                                gwas_locus=pg_tab[['GWAS.locus']][i])
            pg_tab_list[[i]] <- tmp_df
                                    if(i==213) print(tmp_df)
        }
        SNP_gene_pair <- do.call(rbind, pg_tab_list)
    }
    # -
    # Step 2: make grange object from eQTL table
    # -
    if(dataset == 'Science.CD14.Mono'){
        # liftover to hg38, and create a grange object
        grange_SNP <- make_grange_hg19_to_38(SNP_gene_pair, "chr", "snp.pos", "snp.pos")
    }else if (dataset == 'NG.2023.Microglia.GWAS.pcHiC' | dataset == 'NG.2023.Microglia.GWAS.priortized'){
        grange_SNP <- GenomicRanges::makeGRangesFromDataFrame(df = SNP_gene_pair, 
                                                         start.field = 'snp.pos',
                                                         end.field = 'snp.pos',
                                                            keep.extra.columns=TRUE)
    }else if (dataset == 'NG.2023.Microglia.GWAS.pcHiC.peak'){
        grange_SNP <- make_grange_hg19_to_38(SNP_gene_pair, "chr", "peak.start", "peak.end")
    }
    # subset results to only consider the pairs that have genes overlapped with eQTL
    results <- results[results$gene %in% SNP_gene_pair$gene,]
    grange_results <- StringToGRanges(results$peak, sep = c('-', '-'))
    grange_results$gene <- results$gene
    grange_results$pval <- results$pval
    # -
    # Step 3: match SNP and peaks & match genes
    # -
    # number of SNPs in all results pairs (e.g. eQTLs)
    n_pair_all_res <- count_pair_match(grange_results, grange_SNP, 
                                   results$gene, grange_SNP$gene)
    n_pair_all <-  n_pair_all_res[[1]]
    # number of cis-eQTL among significant pairs
    n_pair_sig_res <- count_pair_match(grange_results[results[[pvar]] < p_cutoff], grange_SNP, 
                                   results$gene[results[[pvar]] < p_cutoff], grange_SNP$gene)
    n_pair_sig <- n_pair_sig_res[[1]]
    # number of all pairs that have peaks overlapped with eQTL data
    pval <- phyper(n_pair_sig[2], #n_SNP_sig, 
                   n_pair_all[2], #n_SNP_all, 
                   n_pair_all[1]-n_pair_all[2], #nrow(results) - n_SNP_all, 
                   n_pair_sig[1],#sum(results[[pvar]] < p_cutoff), 
                   lower.tail=F)
    #ratio_sig <- n_SNP_sig / sum(results[[pvar]] < p_cutoff)
    #ratio_bg <- n_SNP_all / nrow(results)
    ratio_sig <- n_pair_sig[2] / n_pair_sig[1]
    ratio_bg <- n_pair_all[2] / n_pair_all[1]
    fc <- ratio_sig / ratio_bg
    return(list(counts = c(n_pair_sig, n_pair_all),
                #c(n_SNP_sig, sum(results[[pvar]] < p_cutoff), n_SNP_all, nrow(results)),
               summary = c(pval, ratio_sig, ratio_bg, fc),
                results[results[[pvar]] < p_cutoff,][n_pair_sig_res[[2]],c('peak','gene','pval')],#,'test_stat'
                n_pair_sig_res[[3]],
               gr_result = grange_results, gr_SNP = grange_SNP, result = results, SNP_gene_pair=SNP_gene_pair,
               matched_peak_gene_gr = grange_results[results[[pvar]] < p_cutoff][n_pair_sig_res[[2]]], 
                matched_SNP_gene_gr = n_pair_sig_res[[3]]))
    
}

TF_motif_enrich <- function(proposed, sparsity, pvar = 'pval', p_cutoff = 0.05, #sparsity_cutoff = 0.5,
                           customized_bg = F){
    TF_not_sparse <- names(which(sparsity[TF_names] < 0.99))
    exp_bg_motifs <- colnames(motif.object@data)[match(TF_not_sparse, full_TF_names)]
    # test enrichment
    # FindMotifs implements a hypergeometric test 
    ct.enriched.motifs <- FindMotifs(
      object = obj,
      features = unique(proposed$peak)#,
      #  background = exp_bg_motifs
    )
    ct.enriched.motifs$sparsity <- sparsity[match(ct.enriched.motifs$motif.name, names(sparsity))]
    #spar1 <- sparsity[ct.enriched.motifs$motif.name] > sparsity_cutoff
    linked.enriched.motifs <- FindMotifs(
      object = obj,
      features = proposed$peak[which(proposed[[pvar]] < p_cutoff)]#,
      #  background = exp_bg_motifs
    )
    linked.enriched.motifs$sparsity <- sparsity[match(linked.enriched.motifs$motif.name, names(sparsity))]
    return(list(ct.enriched.motifs,
                linked.enriched.motifs))
}

make_trio <- function(TF, results, TF_target_p, 
                      TF_peak_p_df, 
                      pair_sigp_cutoff = 0.05, simul = T){
    # extract all peaks that contain the motif for this TF
    peaks_w_binding_motifs <- names(which(motif.object@data[,full_TF_names == TF] != 0))
    results_peaks_w_motif_inds <- results$peak %in% peaks_w_binding_motifs
    if(any(results_peaks_w_motif_inds)){
        # extract two p values, focusing on peaks and genes related to the TF
        # from peak to gene
        peak_gene_pval <- results$pval[results_peaks_w_motif_inds]
        # from gene to TF
        gene_TF_pval <- TF_target_p[TF, results$gene][results_peaks_w_motif_inds]
        # calculate the associations between TF and peak
        peak_TF_pval <- TF_peak_p_df$pval[match(results$peak[results_peaks_w_motif_inds], rownames(TF_peak_p_df))]
        # criterion for a trio
        if(simul){
            trio_inds <- (peak_gene_pval < pair_sigp_cutoff) & (gene_TF_pval < pair_sigp_cutoff) & (peak_TF_pval < pair_sigp_cutoff)
        }else{
            trio_inds <- sapply((peak_gene_pval < pair_sigp_cutoff) + (gene_TF_pval < pair_sigp_cutoff) + (peak_TF_pval < pair_sigp_cutoff),
                               function(x) x >= 2)
        }
        if(any(trio_inds)){
            tmp <- results[results_peaks_w_motif_inds, c('peak','gene','pval')][trio_inds,]
            colnames(tmp)[colnames(tmp) == 'pval'] <- 'pval_peak_gene'
            tmp$pval_TF_gene <- gene_TF_pval[trio_inds]
            tmp$pval_peak_TF <- peak_TF_pval[trio_inds]
            tmp$TF <- TF
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