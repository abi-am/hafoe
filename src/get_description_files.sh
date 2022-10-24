#!/bin/bash

while getopts d:s: flag
do
    case "${flag}" in
        d) needed_wd=${OPTARG};;
        s) samtools_path=${OPTARG};;
    esac
done

if [ ! -d $needed_wd/csv ]; then
  mkdir $needed_wd/csv
fi

if [ ! -d $needed_wd/csv/mapped ]; then
  mkdir $needed_wd/csv/mapped
fi

if [ ! -d $needed_wd/csv/unmapped ]; then
  mkdir $needed_wd/csv/unmapped
fi

echo -ne "\nProcessing mapped reads\n\n"

cnt=0
LEN=$(ls $needed_wd/bam/ | wc -l)
printed0=false
printed25=false
printed50=false
printed75=false
printed100=false
start=$(date +%s)
for f in $needed_wd/bam/*; do
    percent=$((cnt * 100 / LEN))
  
    if [[ $percent -ge 100 ]] && [ "$printed100" = false ] ; then
        #echo -ne "Running...: 100% done\n"
        printed100=true
        end100=$(date +%s)
        echo -ne "\nRunning...: 100% done\nElapsed Time: $(($end100-$start)) seconds\n"
    elif [[ $percent -ge 75 ]] && [ "$printed75" = false ] ; then
        #echo -ne "Running...: 75% done\n"
        printed75=true
        end75=$(date +%s)
        echo -ne "\nRunning...: 75% done\nElapsed Time: $(($end75-$start)) seconds\n"
    elif [[ $percent -ge 50 ]] && [ "$printed50" = false ] ; then
        #echo -ne "Running...: 50% done\n"
        printed50=true
        end50=$(date +%s)
        echo -ne "\nRunning...: 50% done\nElapsed Time: $(($end50-$start)) seconds\n"
    elif [[ $percent -ge 25 ]] && [ "$printed25" = false ] ; then
        #echo -ne "Running...: 25% done\n"
        end25=$(date +%s)
        echo -ne "\nRunning...: 25% done\nElapsed Time: $(($end25-$start)) seconds\n"
        printed25=true
    elif [ "$percent" -ge "0" ] && [ "$printed0" != true ] ; then
        #echo -ne "Running...: 0% done\n"
        end0=$(date +%s)
        echo -ne "\nRunning...: 0% done\nElapsed Time: $(($end0-$start)) seconds\n"
        printed0=true
    fi

    cnt=$((cnt+1))


    f_name=$(basename "$f" .bam)
    $samtools_path view $f | cut -f1,3,12 > $needed_wd/csv/mapped/$f_name.csv; #-f1,3,12 for AS:i:
done
end100=$(date +%s)
echo -ne "\nRunning...: 100% done\nElapsed Time: $(($end100-$start)) seconds\n"



echo -ne "\n\nProcessing unmapped reads\n\n"

cnt=0
LEN=$(ls $needed_wd/bam/ | wc -l)
printed0=false
printed25=false
printed50=false
printed75=false
printed100=false
start=$(date +%s)
for f in $needed_wd/unmapped/*; do
    percent=$((cnt * 100 / LEN))
  
    if [[ $percent -ge 100 ]] && [ "$printed100" = false ] ; then
        #echo -ne "Running...: 100% done\n"
        printed100=true
        end100=$(date +%s)
        echo -ne "\nRunning...: 100% done\nElapsed Time: $(($end100-$start)) seconds\n"
    elif [[ $percent -ge 75 ]] && [ "$printed75" = false ] ; then
        #echo -ne "Running...: 75% done\n"
        printed75=true
        end75=$(date +%s)
        echo -ne "\nRunning...: 75% done\nElapsed Time: $(($end75-$start)) seconds\n"
    elif [[ $percent -ge 50 ]] && [ "$printed50" = false ] ; then
        #echo -ne "Running...: 50% done\n"
        printed50=true
        end50=$(date +%s)
        echo -ne "\nRunning...: 50% done\nElapsed Time: $(($end50-$start)) seconds\n"
    elif [[ $percent -ge 25 ]] && [ "$printed25" = false ] ; then
        #echo -ne "Running...: 25% done\n"
        end25=$(date +%s)
        echo -ne "\nRunning...: 25% done\nElapsed Time: $(($end25-$start)) seconds\n"
        printed25=true
    elif [ "$percent" -ge "0" ] && [ "$printed0" != true ] ; then
        #echo -ne "Running...: 0% done\n"
        end0=$(date +%s)
        echo -ne "\nRunning...: 0% done\nElapsed Time: $(($end0-$start)) seconds\n"
        printed0=true
    fi

    cnt=$((cnt+1))

    f_name=$(basename "$f" .fq)
    grep "@" $f > $needed_wd/csv/unmapped/$f_name.csv;
done
end100=$(date +%s)
echo -ne "\nRunning...: 100% done\nElapsed Time: $(($end100-$start)) seconds\n"
