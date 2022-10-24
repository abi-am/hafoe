#!/bin/bash

# Cluster enriched library sequences with chimeric library representative sequences using cd-hit

while getopts h:i:e:o:p:l:c:n:g:a:m: flag
do
    case "${flag}" in
        h) cd_hit_est_2d_path=${OPTARG};;
        i) input_chimeric_path=${OPTARG};;
        e) input_enriched_path=${OPTARG};;
        o) output_dir=${OPTARG};;
        p) output_prefix=${OPTARG};;
        l) log_path=${OPTARG};;
        c) seq_identity=${OPTARG};;
        n) word_length=${OPTARG};;
        g) identity_type=${OPTARG};;
        a) aln_coverage=${OPTARG};;
        m) mem_limit=${OPTARG};;
    esac
done

echo -ne "Command: $cd_hit_est_2d_path \n\t -i $input_path \n\t -o $output_dir/$output_prefix \n\t -g 1 \n\t -d 100 \n\t -c $seq_identity \n\t -n $word_length \n\t -G $identity_type \n\t -aL $aln_coverage \n" 
echo -ne "cd-hit-est path: $cd_hit_est_2d_path \n"
echo -ne "Input 1 (chimeric representatives) path: $input_chimeric_path \n"
echo -ne "Input 2 (enriched library) path: $input_enriched_path \n"
echo -ne "Output path: $output_dir \n"
echo -ne "Log path: $log_path \n"
echo -ne "Sequence identity threshold, default 0.9 (-c): $seq_identity \n"
echo -ne "Identity type, default 1 (-G):  $identity_type 
\t(if set to 1, use global sequence identity, 
\tif set to 0, then use local sequence identity, 
\tcalculated as number of identical bases in alignment
\tdivided by the length of the alignment) \n"
echo -ne "Alignment coverage for the longer sequence, default 0.0 (-aL): $aln_coverage \n"
echo -ne "word_length, default 10, see CD-HIT user's guide for choosing it (-n): $word_length \n"
echo -ne "memory limit (in MB) for the program, default 800, 0 for unlimitted (-M): $mem_limit \n"


$cd_hit_est_2d_path \
    -i $input_chimeric_path \
    -i2 $input_enriched_path \
    -o $output_dir/$output_prefix \
    -g 1 \
    -d 100 \
    -c $seq_identity \
    -n $word_length \
    -G $identity_type \
    -aL $aln_coverage \
    -M  $mem_limit \
    &>> $log_path

# as an output .fasta and .fasta.clstr files are generated

# To get needed info from the output files, clstr_9595.fasta.clstr text file is used which contains info about which 
#variants are in one cluster. From that file get members_ordered.csv (variants ordered by clusters) and 
#cluster_sizes.csv (sizes of clusters) files using this script:

#get number of rows corresponding to each cluster by using pattern ">C" (each cluster starts with such row: ">Cluster 0")
awk '!/>C/{count++}/>C/&&count{print count; count = 0} END{if (count) print count}' $output_dir/$output_prefix.clstr > $output_dir/cluster_sizes.csv

#get only variant names from the file using sed command to manipulate with strings
cat $output_dir/$output_prefix.clstr | sed -n -e 's/^.*>//p' | sed -n -e 's/\.\.\..*//p' > $output_dir/members_ordered.csv


#get old representatives
cat $output_dir/$output_prefix.clstr | grep "*" | sed -n -e 's/^.*>//p' | sed -n -e 's/\.\.\..*//p' > $output_dir/representatives.csv

