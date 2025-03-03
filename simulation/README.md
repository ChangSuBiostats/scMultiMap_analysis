# Workflow for simulation experiments

1. Run permutation experiments
  - Generate permutation data: `bash generate_data.sh permutation`
  - Run permutation analysis: `bash evaluate_methods.sh permutation`
  - Summarize permutation results: Run `permutation_results_summary.ipynb` for analysis

3. Run simulation experiments (alternative setting)
  - Generate simulation data: `bash generate_data.sh simulation`
  - Run simulation analysis: `bash evaluate_methods.sh simulation`
  - Summarize simulation results: Run `simulation_results_summary.ipynb` for analysis
