#!/bin/bash

# Index
while getopts b:i:o:p:l: flag
do
    case "${flag}" in
        b) bowtie2_build_path=${OPTARG};;
        i) input_path=${OPTARG};;
        o) output_dir=${OPTARG};;
        p) output_prefix=${OPTARG};;
        l) log_path=${OPTARG};;
    esac
done


echo -ne "Command: $bowtie2_build_path \n\t $input_path  \n\t $output_dir/$output_prefix \n"
echo -ne "bowtie2-build path: $bowtie2_build_path \n"
echo -ne "Input path: $input_path \n"
echo -ne "Output path: $output_dir \n"
echo -ne "Log path: $log_path \n"



$bowtie2_build_path $input_path $output_dir/$output_prefix &>> $log_path
