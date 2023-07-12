#!/bin/bash

while getopts b:d:r:s:i:l:L:D: flag
do
    case "${flag}" in
        b) bowtie2_path=${OPTARG};;
        d) needed_wd=${OPTARG};;
        r) ref_index_path=${OPTARG};;
        s) samtools_path=${OPTARG};;
        i) input_dir=${OPTARG};;
        l) log_path=${OPTARG};;
        L) seed_len=${OPTARG};;
        D) extention_num=${OPTARG};;
    esac
done

if [ ! -d $needed_wd/sam ]; then
  mkdir $needed_wd/sam
fi

if [ ! -d $needed_wd/bam ]; then
  mkdir $needed_wd/bam
fi

if [ ! -d $needed_wd/unmapped ]; then
  mkdir $needed_wd/unmapped
fi


echo -ne "Command: $bowtie2_path \n\t -p 4 \n\t -a \n\t --no-unal \n\t -L 30 \n\t -D 2 \n\t -x $ref_index_path "$f" \n\t --un "$needed_wd/unmapped/$f_name.fq" \n\t -S "$needed_wd/sam/$f_name.sam" \n"
echo -ne "bowtie2 path: $bowtie2_path \n"
echo -ne "Reference index path: $ref_index_path \n"
echo -ne "Output path for mapped reads: $needed_wd/bam \n"
echo -ne "Output path for unmapped reads: $needed_wd/unmapped \n"
echo -ne "Log path: $log_path \n"
echo -ne "Seed length (-L): $seed_len \n"
echo -ne "Number of seed extention failed attemps, default 15 (-D): $extention_num 
\t(Up to <int> consecutive seed extension attempts can "fail"
\tbefore Bowtie 2 moves on, using the alignments found so far. 
\tA seed extension "fails" if it does not yield a new best or
\ta new second-best alignment. This limit is automatically 
\tadjusted up when -k or -a are specified.) \n"




cnt=0
LEN=$(ls $needed_wd/$input_dir/ | wc -l)
printed0=false
printed25=false
printed50=false
printed75=false
printed100=false
start=$(date +%s)

for f in $needed_wd/$input_dir/*; do
    #report the progress start
    ##########################
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
    ##########################
    
    
    f_name=$(basename "$f" .fastq)
    echo -ne "\nRunning alignment on $f_name\n" &>> $log_path
    $bowtie2_path \
        -p 4 \
        -a \
        --no-unal \
        -L $seed_len \
        -D $extention_num \
        -x $ref_index_path "$f" \
        --un "$needed_wd/unmapped/$f_name.fq" \
        -S "$needed_wd/sam/$f_name.sam" \
        &>> $log_path
    $samtools_path view -bSh $needed_wd/sam/$f_name.sam > $needed_wd/bam/$f_name.bam
    rm $needed_wd/sam/$f_name.sam
done
rm -r $needed_wd/sam
end100=$(date +%s)
echo -ne "\nRunning...: 100% done\nElapsed Time: $(($end100-$start)) seconds\n"


# cat $needed_wd/log/summary_stat/summary_initial_query.csv | grep -v "Use of uninitialized" > $needed_wd/log/summary_stat/summary_query.csv
# rm $needed_wd/log/summary_stat/summary_initial_query.csv
