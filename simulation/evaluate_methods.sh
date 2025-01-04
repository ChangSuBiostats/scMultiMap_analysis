mode=$1

##
# PBMC data
##
ct="CD14 Mono"
study="PBMC_10k_nextgem"
# scMultiMap
Rscript simulation_proposed.R --ct="$ct" --study=$study --mode=$mode
# Signac
Rscript simulation_Signac.R --ct="$ct" --study=$study --mode=$mode
# SCENT (computationally intensive)
for i in {1..100}
do 
  Rscript permutation_SCENT.R --i_permu=$i --ct="$ct" --study=$study --mode=$mode
done

##
# brain data
##
ct="Oligodendrocytes"
study="brain_CG"
# scMultiMap
Rscript simulation_proposed.R --ct="$ct" --study=$study --mode="$mode"_batch_effect --adjust_for_batch=TRUE
# Signac
Rscript simulation_Signac.R --ct="$ct" --study=$study --mode="$mode"_batch_effect
# SCENT (computationally intensive)
for i in {1..100}
do 
  Rscript permutation_SCENT.R --i_permu=$i --ct="$ct" --study=$study --mode="$mode"_batch_effect --adjust_for_batch=TRUE
done

