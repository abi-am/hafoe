#!/bin/bash
#SBATCH --mem 10gb
#SBATCH --output=slurm-%j.log

#run from hafoe/ directory

./hafoe.sh \
    -parentlib example/input_files/AAV_all16_new.fasta \
    -chimericlib example/input_files/data1/Chimeric_lib_simulated_ccs.csv \
    -enrichedlib1 example/input_files/data1/enriched \
    -o example/hafoe_out \
    -cdhitest ~/cd-hit-v4.8.1-2019-0228/cd-hit-est \
    -cdhitest2d ~/cd-hit-v4.8.1-2019-0228/cd-hit-est-2d \
    -rlib /home/tatevik/R/x86_64-pc-linux-gnu-library/4.1 \
    --explore \
    --identify \
    -title_of_the_run data1_simulated
