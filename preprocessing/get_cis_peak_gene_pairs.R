# Preprocessing script
# Obtain pairs of gene and its nearby peaks to infer associations by different methods

library(dplyr)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)

library("optparse")
option_list = list(
    make_option(c("--study"), type="character", default="brain_CG",
              help="which dataset to analyze", metavar="character"),
    make_option(c("--p_top_gene"), type="integer", default=2000,
              help="top genes to study", metavar="integer"),
    make_option(c("--p_top_peak"), type="integer", default=20000,
              help="top peaks to study", metavar="integer"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
study <- opt$study
p_top_gene <- opt$p_top_gene
p_top_peak <- opt$p_top_peak

print(study)
output_dir <- sprintf('output/%s', study)

# load data
if(study == 'brain_CG'){
    pn <- 'peaks'
    obj_suffix <- '_control'
    cts <- c('Microglia', 'Excitatory', 'Inhibitory', 'Oligodendrocytes', 'Astrocytes')
}else if(grepl('PBMC', study)){
    pn <- 'peaks_ct'
    obj_suffix <- ''
    cts <- c('CD14 Mono', 'CD4 TCM', 'CD8 Naive', 'CD4 Naive', 'CD8 TEM')
}else if(grepl('BMMC', study)){
    pn <- 'peaks'
    obj_suffix <- ''
    cts <- c('CD8+ T', 'CD14+ Mono', 'NK', 'CD4+ T activated', 'Naive CD20+ B', 'Erythroblast', 'CD4+ T naive')
}


# load gene location
source('peak_annotation_helper.R')
gene_500kb_region <- load.promoters(to_left=5e+05, to_right=5e+05)

for(ct in cts){
    pre_dir <- sprintf('output/%s/%s', study, ct)
    # load seurat object 
    ct_obj <- readRDS(sprintf('%s_seurat_obj%s.rds', pre_dir, obj_suffix))
    DefaultAssay(ct_obj) <- pn
    peak_grange <- Signac::granges(x = ct_obj)
    # extract raw counts
    rna_counts <- ct_obj[['RNA']]@counts
    peak_counts <- ct_obj[[pn]]@counts
    # select genes and peaks with high mean levels
    rna_seq_depths <- colSums(rna_counts)
    peak_seq_depths <- colSums(peak_counts)

    gene_mean <- rowSums(rna_counts) / sum(rna_seq_depths)
    peak_mean <- rowSums(peak_counts) / sum(peak_seq_depths)
    gene_order <- order(gene_mean, decreasing = T)
    peak_order <- order(peak_mean, decreasing = T)

    # focus on genes with expression levels ranked top p_top_gene
    top_inds <- gene_order[1:p_top_gene]
    top_genes <- rownames(rna_counts)[top_inds]
    # filter it down to the ones that overlap with annotation
    top_genes <- intersect(top_genes, gene_500kb_region$gene)
    sprintf('percent of genes with annotations: %.2f', length(top_genes) / p_top_gene) %>% print
    # focus on peaks with accessibility levels ranked top p_top_gene
    peak_top_inds <- peak_order[1:p_top_peak]

    ## extract peaks that are close to genes
    # obtain peak gene distance matrix
    peak_distance_matrix <- get_peak_gene_matching(peak_grange, gene_500kb_region)
    peak_distance_matrix <- peak_distance_matrix[, match(top_genes, colnames(peak_distance_matrix))]
    # extract abundant peaks that are near the abundant genes
    peak_candidates <- intersect(which(rowSums(peak_distance_matrix) > 0), peak_top_inds)
    sprintf('precent of abundant peaks with near abundant genes: %.2f', length(peak_candidates) / length(peak_top_inds)) %>% print
    # number of top genes which have abundant peaks nearby
    peak_distance_matrix <- peak_distance_matrix[peak_candidates, ] # subset to only those abundant and nearby peaks
    sprintf('percent of top genes matched: %.2f', mean(colSums(peak_distance_matrix) > 0)) %>% print

    summ <- summary(peak_distance_matrix)
    df <- data.frame(peak = rownames(peak_distance_matrix)[summ$i], 
                     gene = colnames(peak_distance_matrix)[summ$j])
    df$gene_chrom <- seqnames(gene_500kb_region)[match(df$gene, gene_500kb_region$gene)] %>% as.character
    df$gene_mean <- gene_mean[df$gene]
    df$peak_mean <- peak_mean[df$peak]
    write.table(rownames(rna_counts)[top_inds], sprintf('%s/%s_top_%i_genes.txt', output_dir, ct, p_top_gene))
    write.table(rownames(peak_counts)[peak_top_inds], sprintf('%s/%s_top_%i_peaks.txt', output_dir, ct, p_top_peak))
    write.table(df, sprintf('%s/%s_%i_peak_%i_gene_pairs.txt', output_dir, ct, p_top_peak, p_top_gene), quote = F, row.names=F)
    print(sprintf('%s/%s_%i_peak_%i_gene_pairs.txt', output_dir, ct, p_top_peak, p_top_gene))
}
    
    
