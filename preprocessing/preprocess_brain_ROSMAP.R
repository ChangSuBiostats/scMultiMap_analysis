# Preprocessing scripts 
# Analysis: converting insertion counts to fragment counts
# Study: https://pubmed.ncbi.nlm.nih.gov/37774680/
# Data source: https://personal.broadinstitute.org/bjames/AD_snATAC/MFC/
#
#

library(anndata)
library(Seurat)
library(Signac)
library(hdf5r)
library(dplyr)
library(MatrixExtra)
# obtain counts
data_dir <- '/projects/chang/data/ROSMAP/MIT_ROSMAP_Multiomics/bjames/'
ad_rna <- read_h5ad(sprintf("%s/MFC_Multiome_rna_annadata.h5ad", data_dir))
ad_atac <- read_h5ad(sprintf("%s/MFC_Multiome_atac_annadata.h5ad", data_dir))

rna_counts <- ad_rna$layers['raw'] %>% t
peak_counts <- ad_atac$layers['raw'] %>% t
# match cell IDs between RNA and peak data
cell_match_inds <- match(rownames(ad_rna$obs), rownames(ad_atac$obs))
peak_counts <- peak_counts[, cell_match_inds]
all(colnames(peak_counts) == colnames(rna_counts)) %>% print
all(colnames(peak_counts) == rownames(ad_rna$obs)) %>% print

# convert insertion-based counts to fragment counts
# Reference: Section Methods/Fragment computation in https://www.nature.com/articles/s41592-023-02112-6#Sec2
# "The standard 10x protocol for generating the cell-peaks matrix is to count the fragment ends (reads). To estimate fragment counts, we rounded all uneven counts to the next highest even number and halved the resulting read counts."
# https://stackoverflow.com/questions/12649082/r-way-to-map-over-all-entries-in-a-sparse-matrix
summ <- summary(peak_counts)
odd_inds <- which(summ $x %% 2 !=0)
summ$x[odd_inds] <- (summ$x[odd_inds]+1)
summ$x <- summ$x / 2
peak_counts_fragments <- sparseMatrix(i = summ$i, j = summ$j, x = summ$x)
colnames(peak_counts_fragments) <- colnames(peak_counts)
# reformat the peak names
peak_grange <- StringToGRanges(rownames(peak_counts), sep = c(":", "-"))
all_peak_chr <- GRangesToString(grange = peak_grange, sep = c("-", "-"))
rownames(peak_counts_fragments) <- rownames(peak_counts) <- all_peak_chr

ROSMAP <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA"
)
ROSMAP[['peaks']] <- CreateChromatinAssay(peak_counts_fragments)
ROSMAP[['peaks_insertion']] <- CreateChromatinAssay(peak_counts)

ROSMAP@meta.data <- cbind(ad_rna$obs, ad_atac$obs[cell_match_inds,])

for(redu in names(ad_rna$obsm)){
    tmp <- as.matrix(ad_rna$obsm[[redu]])
    rownames(tmp) <- colnames(ROSMAP)
    ROSMAP[[redu]] <- CreateDimReducObject(embeddings = tmp, key = redu, assay = 'RNA')
}

for(redu in names(ad_atac$obsm)){
    tmp <- as.matrix(ad_atac$obsm[[redu]])
    rownames(tmp) <- colnames(ROSMAP)
    if(redu == 'X_uamp') redu <- 'X_umap_ATAC'
    ROSMAP[[redu]] <- CreateDimReducObject(embeddings = tmp, key = redu, assay = 'peaks_insertion')
}

library(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(ROSMAP) <- 'peaks'
# get DNA sequence information for peaks
ROSMAP <- RegionStats(
    object = ROSMAP,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(ROSMAP, 'output/brain_ROSMAP/processed_seurat_obj_brain_ROSMAP.rds')
