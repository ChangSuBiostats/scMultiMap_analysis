load.promoters <- function(gtf="/projects/chang/data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz", 
                           # gtf: A reference file 10x uses to map transcripts to genes
                           # can be downloaded from https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest?
                           # Version: GRCh38 Reference - 2020-A-2.0.0 (May 3, 2021)
                           to_left=2000, to_right=100) {
    # reference: https://github.com/KellisLab/AD_regulome_analysis/blob/main/integration_linking_modules/fit_adgwas.R
    require(GenomicRanges)
    gtf = rtracklayer::readGFF(gtf)
    gtf = gtf[gtf$type == "gene",]
    gtf$tss = gtf$start
    gtf$tss[gtf$strand == "-"] = gtf$end[gtf$strand == "-"]
    gtf$to_left = to_left
    gtf$to_right = to_right
    gtf$to_left[gtf$strand == "-"] = to_right
    gtf$to_right[gtf$strand == "-"] = to_left
    d = duplicated(gtf$gene_name)
    gtf$gene_name[d] = paste0(gtf$gene_name[d], "-1")
    gr = with(gtf, GRanges(seqnames=seqid,
                           ranges=IRanges(tss - to_left, tss + to_right),
                           gene=gene_name,
                           tss=tss,
                           strand=strand))
    names(gr) = gr$gene
    return(gr)
}

get_peak_gene_matching <- function(peaks, genes, sep = c("-", "-")){
    # reference: https://rdrr.io/cran/Signac/src/R/links.R
    overlaps <- GenomicRanges::findOverlaps(query = peaks, subject = genes,
        type = "any", select = "all")
    hit_matrix <- Matrix::sparseMatrix(i = S4Vectors::queryHits(x = overlaps),
        j = S4Vectors::subjectHits(x = overlaps), x = 1, dims = c(length(x = peaks),
            length(x = genes)))
    rownames(x = hit_matrix) <- GRangesToString(grange = peaks,
        sep = sep)
    colnames(x = hit_matrix) <- genes$gene
    return(hit_matrix)
}

get_grange_from_rn <- function(rn){
    rn_c <- character(length(rn))
    for(i in 1:length(rn)){
        tmp <- strsplit(rn[i], '\\-')[[1]]
        rn_c[i] <- paste(paste(tmp[1], tmp[2], sep = ':'), tmp[3], sep='-')
    }
    GRanges(rn_c)
}