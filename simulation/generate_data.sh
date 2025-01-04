mode=$1

# PBMC data
ct="CD14 Mono"
study="PBMC_10k_nextgem"
if [ "$mode" = "simulation" ]
then  
  Rscript simulation_realistic_generate_data.R --ct="$ct" --study="$study" --simulation_mode=simulation
else
  Rscript permutation_generate_data.R --ct=$ct --study="$study"
fi

# brain data
ct="Oligodendrocytes"
study="brain_CG"
if [ "$mode" = "simulation" ]
then  
  Rscript simulation_realistic_generate_data.R --ct="$ct" --study="$study" --simulation_mode=simulation_batch_effect
else
  Rscript permutation_generate_data.R --ct=$ct --study="$study" --permu_within_batch=TRUE
fi
