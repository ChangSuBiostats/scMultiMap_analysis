# subset seurat object to cell types
# and save the subsetted object
library(Seurat)
library("optparse")
option_list = list(
    make_option(c("--study"), type="character", default="brain_CG",
              help="which dataset to analyze", metavar="character"),
    make_option(c("--cell_subset"), type="character", default="Control",
		help="which subset of cells to analyze", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
study <- opt$study
cell_subset <- opt$cell_subset

print(study)
print(cell_subset)

# load data
if(study == 'brain_CG'){
    data_dir <- '/projects/chang/scGRN/peak_gene'
    obj <- readRDS(sprintf('%s/data/brain/processed_AD_DLPFC_15.rds', data_dir))
    # subset by disease status and cell type
    if(cell_subset == 'Control'){
        control_inds <- (obj$Diagnosis == 'Unaffected')
        obj_suffix <- '_control'
    }else if(cell_subset == 'All'){
        control_inds <- rep(T, nrow(obj))
        obj_suffix <- ''
    }
    ct_var <- obj$predicted.id
    ct_var_name <- 'predicted.id'
    pn <- 'peaks'
    cts <- c('Microglia', 'Excitatory', 'Inhibitory', 'Oligodendrocytes', 'Astrocytes')
}else if(grepl('PBMC', study)){
    # load data preprocessed with preprocess_10x_PBMC_data.R
    data_dir <- '/projects/chang/scGRN/data/PBMC'
    obj <- readRDS(sprintf('output/PBMC/processed_seurat_obj_%s.rds', study))
    # obj <- readRDS('/projects/chang/scGRN/signac/PBMC/output/PBMC_seurat_obj_peaks_by_ct.rds')
    # obj <- readRDS('/projects/chang/scGRN/signac/PBMC/output/PBMC_seurat_obj_processed.rds')
    control_inds <- rep(T, ncol(obj))
    ct_var <- obj$predicted.id
    ct_var_name <- 'predicted.id'
    pn <- 'peaks_ct'
    cts <- c('CD14 Mono', 'CD4 TCM', 'CD8 Naive', 'CD4 Naive', 'CD8 TEM')
}else if(grepl('BMMC', study)){
    data_dir <- '/projects/chang/scGRN/data/BMMC'
    obj <- readRDS(sprintf('output/BMMC/processed_seurat_obj_%s.rds', study))
    control_inds <- rep(T, ncol(obj))
    ct_var <- obj$cell_type
    pn <- 'peaks'
    cts <- c('CD8+ T', 'CD14+ Mono', 'NK', 'CD4+ T activated', 'Naive CD20+ B', 'Erythroblast', 'CD4+ T naive')
}else if(study == 'brain_ROSMAP'){
    obj <- readRDS('output/brain_ROSMAP/processed_seurat_obj_brain_ROSMAP.rds')
    # subset by disease status and cell type
    if(cell_subset == 'Control'){
        control_inds <- (obj$Pathology == 'nonAD')
        obj_suffix <- '_control'
    }else if(cell_subset == 'All'){
        control_inds <- rep(T, nrow(obj))
        obj_suffix <- ''
    }
    ct_var <- obj$MajorCellType
    ct_var_name <- 'MajorCellType'
    pn <- 'peaks'
    cts <- c('Mic', 'Exc', 'Inh', 'Oli', 'Ast')
}

for(ct in cts){
    print(ct)
    if(study == 'brain_CG'){
        ct_obj <- subset(x = obj, subset = (predicted.id == ct) & control_inds)
    }else if(grepl('PBMC', study)){
        ct_obj <- subset(x = obj, subset = (predicted.id == ct))
        obj_suffix <- ''
    }else if(grepl('BMMC', study)){
        ct_obj <- subset(x = obj, subset = (cell_type == ct))
        obj_suffix <- ''
    }else if(study == 'brain_ROSMAP'){
        ct_obj <- subset(x = obj, subset = (MajorCellType == ct) & control_inds)
    }

    # evaluate SCTransform residuals to compute Signac correlations
    set.seed(2024)
    ct_obj <- Seurat::SCTransform(ct_obj, vst.flavor="v2")

    if(study == 'brain_CG'){
        library(BSgenome.Hsapiens.UCSC.hg38)
    DefaultAssay(ct_obj) <- pn
        # get DNA sequence information for peaks
        ct_obj <- RegionStats(
            object = ct_obj,
            genome = BSgenome.Hsapiens.UCSC.hg38,
            sep = c(":", "-")
        )
    }

    saveRDS(ct_obj, 
            sprintf('output/%s/%s_seurat_obj%s.rds', study, ct, obj_suffix))
    print(sprintf('output/%s/%s_seurat_obj%s.rds', study, ct, obj_suffix))
}
    
