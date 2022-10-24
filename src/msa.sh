#!/bin/bash

# Multiple sequence alignment using cd-hit

while getopts c:i:o: flag
do
    case "${flag}" in
        c) clustalo_path=${OPTARG};;
        i) input_path=${OPTARG};;
        o) output_path=${OPTARG};;

    esac
done

echo -ne "Command: $clustalo_path \n\t -i $input_path \n\t -o $output_path \n\t --outfmt clu \n" 
echo -ne "clustalo path: $clustalo_path \n"
echo -ne "Input path: $input_path \n"
echo -ne "Output path: $output_path \n"


$clustalo_path -i $input_path -o $output_path --outfmt clu --force
