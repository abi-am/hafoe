#!/bin/bash

# Run hafoe on 5 datasets with all combinations of input parameters having following values:
# -readlength - 50,100,150,200
# -stepsize - 10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200
# -vd_criterion - "avg", "sum"

# Create the output folder
mkdir benchmarking 

for d in {1..5}
do
    for rl in {50..200..50}
    do
	for ss in {10..200..10}
	do
	    if [ "$rl" -ge "$ss" ]
	    then
		declare -a arr=("avg" "sum")
		for criterion in "${arr[@]}"
		do
		    ./hafoe.sh \
				-parentlib input_files/AAV_all16_new.fasta \
				-chimericlib input_files/data${d}/Chimeric_lib_simulated_ccs.csv \
				-enrichedlib1 input_files/data${d}/Enriched_lib_simulated_ccs.fastq.gz  \
				-o benchmarking/hafoe_out_${d}_${rl}_${ss}_${criterion} \
				-cdhitest ~/cd-hit-v4.8.1-2019-0228/cd-hit-est \
				-cdhitest2d ~/cd-hit-v4.8.1-2019-0228/cd-hit-est-2d \
				-rlib /home/tatevik/R/x86_64-pc-linux-gnu-library/4.1 \
				--explore \
				-readlength $rl \
				-stepsize $ss \
				-vd_criterion $criterion \ 
				-title_of_the_run simulated_data${d}_${rl}_${ss}_${criterion}
		done
	    fi
	done
    done
done