library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

library("optparse")
option_list = list(
    make_option(c("--study"), type="character", default="brain_CG",
              help="which dataset to analyze", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
study <- opt$study

# load full single cell object for all peaks
data_dir <- '/projects/chang/scGRN/peak_gene'
if(study == 'brain_CG'){
    obj <- readRDS(sprintf('%s/data/brain/processed_AD_DLPFC_15.rds', data_dir))
    pn <- "peaks"
}else if(study == 'PBMC_4_combined'){
    obj <- readRDS(sprintf('../preprocessing/output/PBMC/processed_seurat_obj_%s.rds', 'PBMC_4_combined'))
    pn <- 'peaks_ct'
}
motif_fn <- sprintf('output/%s/motif.object.rds', study)
DefaultAssay(obj) <- pn

if(!file.exists(motif_fn)){
    print(sprintf('Generate motif object for study %s', study))
    # Get a list of motif position frequency matrices from the JASPAR database
    pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
    # add motif information
    motif.matrix <- CreateMotifMatrix(features = granges(obj), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
    motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
    saveRDS(motif.object, motif_fn)
}