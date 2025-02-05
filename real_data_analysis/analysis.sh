# Run all three methods for reproducibility, consistency, and trio analysis.
# Note that
# 1. For Signac, when multiple subjects are present in the data, we further adjusted Signac by its estimates on permuted null data that preserve across-subject variations. Please refer to subsection "Reproducibility analysis" in the manuscript for more details.
# 2. For SCENT, we split all peak-gene pairs into multiple sets due to extreme computational cost.

###############################################################
# 1. PBMC
###############################################################

# -------------------------------------------------------------
# 1.1 reproducibility across technical and biological replicates
# -------------------------------------------------------------
for study in PBMC_10k_nextgem PBMC_10k_nextgem_X PBMC_10k
do
    # run scMultiMap
    Rscript run_proposed.R --i_ct=1 --study=$study 
    # run Signac
    Rscript run_Signac.R --i_ct=1 --study=$study 
    # run SCENT; split the peak-gene pairs into 10 sets due to extreme computational cost
    n_set=10
    for i_set in {1..10}
    do
        Rscript run_SCENT.R --ct='CD14 Mono' --study=$study --n_set=10 --i_set=$i_set  --n_bootstrap=5000
    done
done

# -------------------------------------------------------------
# 1.2 consistency with orthogonal data modalities
# -------------------------------------------------------------

# run scMultiMap
Rscript run_proposed.R --i_ct=1 --study=PBMC_4_combined --adjust_for_batch=TRUE

# run Signac
Rscript run_Signac.R --i_ct=1 --study=PBMC_4_combined
## adjust Signac with estimates on permuted data
Rscript run_Signac.R --i_ct=1 --study=PBMC_4_combined 
for i in {1..100}
do
    Rscript run_Signac.R --i_ct=1 --study=PBMC_4_combined --permu_within_batch=TRUE --i_permu=$i 
done
Rscript adjust_Signac.R --i_ct=1 --study=PBMC_4_combined

# run SCENT; split the peak-gene pairs into 10 sets due to extreme computational cost
n_set=10
for i_set in {1..10}
do
    Rscript run_SCENT.R --ct='CD14 Mono' --study=$study --adjust_for_batch=TRUE --n_set=10 --i_set=$i_set  --n_bootstrap=5000
done

# -------------------------------------------------------------
# 1.3 trio analysis
# -------------------------------------------------------------

# run scMultiMap
Rscript run_scMultiMap.R --i_ct=1 --study=PBMC_4_combined --adjust_for_batch=TRUE --run_TF_gene=TRUE
Rscript construct_trios.R --i_ct=1 --study=PBMC_4_combined --adjust_for_batch=TRUE 

###############################################################
# 2. Brain
###############################################################

# -------------------------------------------------------------
# 2.1 reproducibility across two studies & consistency with orthogonal data modalities
# -------------------------------------------------------------
for study in brain_CG brain_ROSMAP
do
    if [[ "$study" == "brain_CG" ]]; then
        cts=('Oligodendrocytes' 'Excitatory' 'Inhibitory' 'Microglia' 'Astrocytes')
    else
        cts=('Oli' 'Exc' 'Inh' 'Mic' 'Ast')
    fi

    for i_ct in {1..5}
    do
        # run scMultiMap
        Rscript run_scMultiMap.R --i_ct=$i_ct --study=$study --adjust_for_batch=TRUE --cell_subset=control
        
        # run Signac; adjust Signac with estimates on permuted data
        Rscript run_Signac.R --i_ct=$i_ct --study=$study
        for i in {1..100}
        do
            Rscript run_Signac.R --i_ct=$i_ct --study=$study --permu_within_batch=TRUE --i_permu=$i 
        done
        Rscript adjust_Signac.R --i_ct=$i_ct --study=$study

        # run SCENT; split the peak-gene pairs into 10 sets due to extreme computational cost
        for i_set in {1..20}
        do
            Rscript run_SCENT.R --ct="${cts[${i_ct}]}" --study=$study --n_set=20 --i_set=$i_set --adjust_for_batch=TRUE --n_bootstrap=5000
        done
    done
done

# -------------------------------------------------------------
# 2.2 trio analysis
# -------------------------------------------------------------
for i_ct in {1..5}
do
    # run scMultiMap
    Rscript run_scMultiMap.R --i_ct=$i_ct --study=brain_CG --adjust_for_batch=TRUE --run_TF_gene=TRUE
    Rscript construct_trios.R --i_ct=$i_ct --study=brain_CG --adjust_for_batch=TRUE 
done

# -------------------------------------------------------------
# 2.3 differential co-expression
# -------------------------------------------------------------

for i_ct in 1 5
do
    # obtain scMultiMap estimates using cells from AD subjects
    Rscript run_scMultiMap.R --i_ct=$i_ct --study=brain_CG --adjust_for_batch=TRUE --run_disease=TRUE
    # obtain permutation p-values
    for i_permu in {1..100}
    do
        Rscript run_proposed_differential.R --i_ct=$i_ct --study=brain_CG --n_permu=100 --i_permu=$i_permu
    done
done


