# Workflow for simulation experiments


1. Run permutation experiments

  1. Generate permutation data: `bash generate_data.sh permutation`

  2. Run permutation analysis: `bash evaluate_methods.sh permutation`

  3. Summarize permutation results: Run `permutation_results_summary.ipynb` for analysis

2. Run simulation experiments (alternative setting)

  1. Generate simulation data: `bash generate_data.sh simulation`

  2. Run simulation analysis: `bash evaluate_methods.sh simulation`

  3. Summarize simulation results: Run `simulation_results_summary.ipynb` for analysis
