# Preprocessing scripts 
# Analysis: peak calling by cell type with MACS2
# Data: single cell multiome data on brain from an AD study: https://doi.org/10.1016/j.xgen.2023.100263
# Data source: GSE214637
#
# Code reference: Signac tutorial: 
# [Calling peaks](https://stuartlab.org/signac/articles/peak_calling). 


library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)


# load the RNA and ATAC data
data_dir <- '/projects/chang/data/multiome/AD_DLPFC_10x'
counts <- Read10X_h5(sprintf("%s/GSE214979_filtered_feature_bc_matrix.h5", data_dir))
fragpath <- sprintf("%s/GSE214979_atac_fragments.tsv.gz", data_dir)
metadata <- read.csv(sprintf('%s/GSE214979_cell_metadata.csv', data_dir), row.names=1)

# set directory to save output
output_dir <- '/projects/chang/scGRN/peak_gene'

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
AD <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
AD[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

AD <- AddMetaData(AD, metadata, col.name = NULL)


# -
# call peaks by cell type with MACS2
# -

DefaultAssay(AD) <- "ATAC"
# call peaks using MACS2 by cell type
peaks <- CallPeaks(AD, group.by = 'predicted.id')
warnings() %>% print

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(AD),
  features = peaks,
  cells = colnames(AD)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
AD[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# get DNA sequence information
AD[['peaks']] <- RegionStats(
  object = AD[['peaks']],
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

saveRDS(AD, sprintf('%s/data/brain/processed_AD_DLPFC_15.rds', output_dir))
