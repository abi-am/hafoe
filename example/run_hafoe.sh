#!/bin/bash

# Run from hafoe/ directory

./hafoe.sh \
    -parentlib example/input_files/AAV_all16_new.fasta \
    -chimericlib example/input_files/data1/Chimeric_lib_simulated_ccs.csv \
    -enrichedlib1 example/input_files/data1/enriched \
    -o example/hafoe_out \
    --explore \
    --identify \
    -title_of_the_run data1_simulated
