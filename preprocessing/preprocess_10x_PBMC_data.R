# Preprocessing scripts 
# Analysis: cell type annotation + peak calling by cell type with MACS2
# Data: PBMC multiome data from 10x Genomics. Scripts for downloading data are included below.
#
#
# Code reference: two Signac tutorials: 
# 1. [Joint RNA and ATAC analysis: 10x multiomic](https://stuartlab.org/signac/articles/pbmc_multiomic.html) 
# 2. [Calling peaks](https://stuartlab.org/signac/articles/peak_calling). 
# The first tutorial is for annotating PBMC cell types, and the second tutorial is for calling peaks by cell type to identify as many cell-type-specific peaks as possible. 
#
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

library("optparse")
option_list = list(
    make_option(c("--dataset"), type="character", default="PBMC_3k",
              help="which dataset to analyze", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dataset <- opt$dataset
print(dataset)

data_dir <- '/projects/chang/scGRN/data/PBMC'

# here * represents filtered_feature_bc_matrix.h5/atac_fragments.tsv.gz/atac_fragments.tsv.gz.tbi
dataset_str <- c(
    # wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_*
    'pbmc_granulocyte_sorted_3k', 
    # wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_*
    'pbmc_granulocyte_sorted_10k', 
    # wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_Controller/10k_PBMC_Multiome_nextgem_Chromium_Controller_*
    '10k_PBMC_Multiome_nextgem_Chromium_Controller',
    # wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X*
    '10k_PBMC_Multiome_nextgem_Chromium_X')

names(dataset_str) <- c('PBMC_3k', 'PBMC_10k', 'PBMC_10k_nextgem', 'PBMC_10k_nextgem_X')

# load the RNA and ATAC data
counts <- Read10X_h5(sprintf("%s/%s_filtered_feature_bc_matrix.h5", data_dir, dataset_str[dataset]))
fragpath <- sprintf("%s/%s_atac_fragments.tsv.gz", data_dir, dataset_str[dataset])

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

print(pbmc)

# -
# Quality control
# -
DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

g <- VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
ggsave(sprintf('figures/PBMC/QC_%s.png', dataset), g)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
print('\nPost QC\n')
print(pbmc)

# -
# Peak calling
# -
print('Call peaks')
print(Sys.time())
# call peaks using MACS2
peaks <- CallPeaks(pbmc)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

print('Peak calling ended')
print(Sys.time())

# -
# Gene expression data processing
# -
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

# -
# DNA accessibility data processing
# -
DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

# -
# Annotating cell types
# -
# "We’ll use an annotated PBMC reference dataset from Hao et al. (2020), available for download here: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"

library(SeuratDisk)

# load PBMC reference
reference <- LoadH5Seurat(sprintf("%s/pbmc_multimodal.h5seurat", data_dir), assays = list("SCT" = "counts"), reductions = 'spca')
reference <- UpdateSeuratObject(reference)

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

saveRDS(pbmc, sprintf('output/PBMC/intermediate_seurat_%s.rds', dataset))

print(levels(pbmc))
# set a reasonable order for cell types to be displayed when plotting
full_levels <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                 "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                 "NK Proliferating", "gdT",
                 "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                 "CD14 Mono", "CD16 Mono",
                 "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")
levels(pbmc) <- full_levels[full_levels %in% levels(pbmc)]
print(levels(pbmc))

# -
# Call peaks by cell type
# -
print('Call peaks by cell type')
print(Sys.time())

DefaultAssay(pbmc) <- "ATAC"
# call peaks using MACS2 by cell type
peaks <- CallPeaks(pbmc, group.by = 'predicted.id')
warnings() %>% print

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks_ct"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# get DNA sequence information
pbmc[['peaks_ct']] <- RegionStats(
  object = pbmc[['peaks_ct']],
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

print('peak calling by cell types ended')
print(Sys.time())

saveRDS(pbmc, sprintf('output/PBMC/processed_seurat_obj_%s.rds', dataset))
