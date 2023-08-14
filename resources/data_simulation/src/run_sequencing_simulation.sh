#!/bin/bash

####################################################################
# Simulate PacBio ccs sequencing by simlord (SimLoRD v1.0.4) tool
# using generated counts as number of reads to be generated
####################################################################
# Note: for running this script a miniconda3 env is required,
#       follow the instructions in simlord tool's documentation (https://bitbucket.org/genomeinformatics/simlord/src/master/)
#       and run this script accordingly


for d in {1..5}
do
    echo "Running on data${d}: Chimeric"
    simulate_pacbio_sequencing.sh \
	-i input_files/data${d}/Chimeric_lib_simulated_initial.fasta \
	-c files/data${d}/chimeric_counts.csv \
	-o input_files/data${d}/Chimeric_lib_simulated_ccs.fastq 
done

for d in {1..5}
do
    echo "Running on data${d}: Enriched"
    simulate_pacbio_sequencing.sh \
	-i input_files/data${d}/Chimeric_lib_simulated_initial.fasta \
	-c files/data${d}/enriched_counts.csv \
	-o input_files/data${d}/Enriched_lib_simulated_ccs.fastq
done