# We copied the following codes from https://github.com/stuart-lab/signac/blob/8ecdde2/R/links.R line 225-485
# to implement Signac LinkPeaks()
# while we made three customizations:
# 1. We added the mean and sd of the randomly sampled background distribution as an output of the function
# 2. We fixed the bug in evaluating two-sided p-values (https://github.com/stuart-lab/signac/blob/8ecdde2/R/links.R#L481C29-L481C30)
# 3. We added set.seed() to ensure reproducibility of randomly sampled background distribution
#
# We show that the score and pvalues from this cutomiazed implementation and the original codes by Signac R package
# are the same at sanity_check_Signac_implementations.R
#

# Calculate the correlations between a gene and its neighbouring peaks 'peaks.test'
# and a background distribution that pairs the gene with randomly sampled peaks from other chromosome
Signac_LinkPeaks <- function(
    peaks.test, # peaks to be tested for a gene
    gene.chrom, # which chromosome the selected gene is on, format: chr10
    gene.expression, # normalized gene expression from sctransform slot data
    object, # seurat object for a cell type
    peak.assay = 'peaks', # from the usage in https://stuartlab.org/signac/articles/pbmc_multiomic
    peak.slot = 'counts', # from the default in https://stuartlab.org/signac/reference/linkpeaks
    features.match = c("GC.percent", "count", "sequence.length"), # https://github.com/stuart-lab/signac/blob/8ecdde2/R/links.R#L260
    n_sample = 200,
    seed = 1,
    verbose=F,
    min.cell = 10
){
    peak.data <- GetAssayData(
        object = object, assay = peak.assay, layer = peak.slot
    )
    
    # background parameters
    all.peaks <- rownames(peak.data)
    all.peaks <- all.peaks[rowSums(peak.data > 0) > min.cell] # remove peaks with extremely low counts
    
    meta.features <- GetAssayData(
        object = object, assay = peak.assay, layer = "meta.features"
    )
    
    # select peaks at random with matching GC content, accessibility and peak length
    # sample from peaks on a different chromosome to the gene
    #peaks.test <- rownames(x = coef.result)
    trans.peaks <- all.peaks[
        !grepl(pattern = paste0("^", gene.chrom), x = all.peaks)
    ]
    meta.use <- meta.features[trans.peaks, ]
    pk.use <- meta.features[peaks.test, , drop=FALSE]
    # extract peaks with matched features
    set.seed(seed) ## added
    bg.peaks <- lapply(
            X = seq_len(length.out = nrow(x = pk.use)),
            FUN = function(x) {
              MatchRegionStats(
                meta.feature = meta.use,
                query.feature = pk.use[x, , drop = FALSE],
                features.match = features.match,
                n = n_sample,
                verbose = FALSE
              )
            }
          )
    
    # run background correlations
    peak.data <- t(x = peak.data)
    bg.access <- peak.data[, unlist(x = bg.peaks), drop = FALSE]
    bg.coef <- qlcMatrix::corSparse(
        X = bg.access,
        Y = gene.expression
    )
    

    bg.sum.mat <- sapply(1:length(peaks.test), function(i){
        inds <- ((i-1)*200+1) : (i*200)
        return(c(mean(bg.coef[inds,1]), sd(bg.coef[inds,1])))
    }) %>% t

    # correlations of observed data
    peak.access <- peak.data[, peaks.test, drop = FALSE]
    coef.result <- qlcMatrix::corSparse(
          X = peak.access,
          Y = gene.expression
        )

    df <- data.frame(bg.sum.mat)
    colnames(df) <- c('bg.mean', 'bg.sd')
    df$score <- coef.result
    df$peaks <- peaks.test
    return(df)
}

Signac_LinkPeaks_get_bg_peaks <- function(
    peaks.test, # peaks to be tested for a gene
    gene.chrom, # which chromosome the selected gene is on, format: chr10
    gene.expression, # normalized gene expression from sctransform slot data
    object, # seurat object for a cell type
    peak.assay = 'peaks', # from the usage in https://stuartlab.org/signac/articles/pbmc_multiomic
    peak.slot = 'counts', # from the default in https://stuartlab.org/signac/reference/linkpeaks
    features.match = c("GC.percent", "count", "sequence.length"), # https://github.com/stuart-lab/signac/blob/8ecdde2/R/links.R#L260
    n_sample = 200,
    seed = 1,
    verbose=F,
    min.cell = 10
){
    peak.data <- GetAssayData(
        object = object, assay = peak.assay, layer = peak.slot
    )
    
    # background parameters
    all.peaks <- rownames(peak.data)
    all.peaks <- all.peaks[rowSums(peak.data > 0) > min.cell] # remove peaks with extremely low counts
    
    meta.features <- GetAssayData(
        object = object, assay = peak.assay, layer = "meta.features"
    )
    
    # select peaks at random with matching GC content, accessibility and peak length
    # sample from peaks on a different chromosome to the gene
    #peaks.test <- rownames(x = coef.result)
    trans.peaks <- all.peaks[
        !grepl(pattern = paste0("^", gene.chrom), x = all.peaks)
    ]
    meta.use <- meta.features[trans.peaks, ]
    pk.use <- meta.features[peaks.test, , drop=FALSE]
    # extract peaks with matched features
    set.seed(seed) ## added
    bg.peaks <- lapply(
            X = seq_len(length.out = nrow(x = pk.use)),
            FUN = function(x) {
              MatchRegionStats(
                meta.feature = meta.use,
                query.feature = pk.use[x, , drop = FALSE],
                features.match = features.match,
                n = n_sample,
                verbose = FALSE
              )
            }
          )
    # run background correlations
    peak.data <- t(x = peak.data)
    bg.access <- peak.data[, unlist(x = bg.peaks), drop = FALSE]
    bg.coef <- qlcMatrix::corSparse(
        X = bg.access,
        Y = gene.expression
    )
    return(list(bg.peaks, bg.coef))
}