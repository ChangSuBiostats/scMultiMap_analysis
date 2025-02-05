# follow
# 1. https://stuartlab.org/signac/articles/merging
#
# to merge 4 single cell multiome datasets
library(dplyr)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
datasets <- c('PBMC_3k', 'PBMC_10k', 'PBMC_10k_nextgem', 'PBMC_10k_nextgem_X')

multiome_obj <- list()
for(dataset in datasets){
    multiome_obj[[dataset]] <- readRDS(sprintf('output/PBMC/processed_seurat_obj_%s.rds', dataset))
}

# -
# merge peak data
# -

# Step 1: create a common peak set
grange_list <- list()
for(dataset in datasets){
    grange_list[[dataset]] <- Signac::granges(multiome_obj[[dataset]][['peaks_ct']])
}

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(grange_list[[1]], grange_list[[2]], grange_list[[3]], grange_list[[4]]))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
summary(peakwidths) # all between 200 and 4088
#combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
#combined.peaks

# Step 2: create fragment objects
# skipped as it was part of the original datasets provided by 10x

# Step 3: Quantify peaks in each dataset
common_peak_counts <- list()
for(dataset in datasets[2:4]){
    common_peak_counts[[dataset]] <- FeatureMatrix(
        fragments = Fragments(multiome_obj[[dataset]]),
        features = combined.peaks,
        cells = colnames(multiome_obj[[dataset]])
    )
}

# Step 4: Create new Seurat objects
new_multiome_obj <- list()
for(dataset in datasets){
    # customize: initiate the new seurat object with RNA counts
    new_multiome_obj[[dataset]] <- CreateSeuratObject(
        counts = multiome_obj[[dataset]]$RNA@counts,
        assay = "RNA"
    )
    ChromAssay <- CreateChromatinAssay(common_peak_counts[[dataset]], 
                                       fragments = Fragments(multiome_obj[[dataset]]))
    new_multiome_obj[[dataset]][['peaks_ct']] <- ChromAssay
    new_multiome_obj[[dataset]]$predicted.id <- multiome_obj[[dataset]]$predicted.id
    new_multiome_obj[[dataset]]$dataset <- dataset
}

# Step 5: Merge objects
# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = new_multiome_obj[[1]],
  y = list(new_multiome_obj[[2]], new_multiome_obj[[3]], new_multiome_obj[[4]]),
  add.cell.ids = datasets
)
combined[["peaks_ct"]]

# get DNA sequence information
DefaultAssay(combined) <- 'peaks_ct'

combined <- RegionStats(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

# calculate dimension reduction for peak data
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

# calculate dimension reduction for RNA data
DefaultAssay(combined) <- "RNA"
pbmc <- SCTransform(combined)
pbmc <- RunPCA(combined)

saveRDS(combined, sprintf('output/PBMC/processed_seurat_obj_%s.rds', 'PBMC_4_combined'))