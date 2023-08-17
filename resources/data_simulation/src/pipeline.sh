#!/bin/bash

# Run this script from data_simulation directory

# 1. Run data_simulation_before.R to generate sequences and counts for 5 datasets
# Outputs: input_files/dataX/Chimeric_lib_simulated_initial.fasta
#          input_files/dataX/Chimeric_lib_simulated_labels.csv
#          files/dataX/chimeric_counts.csv
#          files/dataX/enriched_counts.csv
echo "######## RUNNING data_simulation_before.R #######"
Rscript src/data_simulation_before.R

# 2. Run simulate_sequencing_pacbio.sh to get PacBio ccs reads of the generated sequences with corresponding counts
# Outputs: input_files/dataX/Chimeric_lib_simulated_ccs.fastq 
#          input_files/dataX/Enriched_lib_simulated_ccs.fastq            
echo "######## RUNNING run_sequencing_simulation.sh #######"
bash src/run_sequencing_simulation.sh

# 3. Run data_simulation_after.R to generate chimeric csv file after simlord simulation
# Outputs: input_files/dataX/Chimeric_lib_simulated_ccs.csv
echo "######## RUNNING data_simulation_after.R #######"
Rscript src/data_simulation_after.R

# 4. Run hafoe on 5 simulated datasets with different combinations of input parameters
# Outputs: benchmarking/hafoe_out_${d}_${rl}_${ss}_${criterion}
echo "######## RUNNING run_hafoe_benchmarking.sh #######"
bash src/run_hafoe_benchmarking.sh

# 5. Make accuracy comparison boxplot on hafoe benchmarking results
# Outputs: plots/variant_composition_accuracy_comparison.pdf
echo "######## RUNNING accuracy_comparison_plot.R #######"
Rscript src/accuracy_comparison_plot.R

# 6. Multiple seqeunce alignment of top20 representatives for true/predicted composition plot
echo "######## RUNNING get_files_for_true_pred_composition_plot.R #######"
Rscript src/get_files_for_true_pred_composition_plot.R

# 7. Plot true vs predicted composition and positional abundance plots
# see src/true_variant_composition_plot.ipynb notebook


