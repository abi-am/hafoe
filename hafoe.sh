#!/bin/bash

usage="\nProgram:\thafoe
\nVersion:\t1.0.0
\n\nusage:\t./hafoe.sh [options] -parentlib <path_to_parental_sequences_fasta_file> -chimericlib <path_to_chimeric_sequences_csv/fastq_file> -enrichedlib1 <path_to_dir_with_enriched_sequences_fastq_file(s)> -o <output_directory>
\n\nInput:
\n\t-parentlib\tfasta file (required when --explore specified: the full path to parental sequences)
\n\t-chimericlib\tcsv or fastq file (required when --explore specified: the full path to chimeric sequences)
\n\t-enrichedlib1\tdirectory with fastq file(s) (required when --identify specified: the full path to the directory with sequences after 1st enrichment cycle of one or more samples  e.g. your_fq_folder)
\n\t-enrichedlib2\tdirectory with fastq file(s) (optional: the full path to the directory with sequences after 2nd enrichment cycle,  --identify should be specified)
\n\t-exploreout\t (required when only --identify specified: the full path to output directory generated by running hafoe with only --explore option)

\n\nOptions:
\n\t--explore\toption directs hafoe to run the first part of the program: exploring the chimeric dataset
\n\t--identify\toption directs hafoe to run the second part of the program: identifying successful variants after enrichment cycles
\n\tAt least one of --explore, --identify options should be specified.
\n\t-o\toutput directory (optional: the default is hafoe_out)
\n\t-samtools\tsamtools path (optional: if not supplied, hafoe will use the samtools installed on the system)
\n\t-bowal\tbowtie2 path (optional: hafoe will use the default installation path in the user's directory)
\n\t-bowb\tbowtie2-build path (optional: hafoe will use the default installation path in the user's directory)
\n\t-cdhitest\tcd-hit-est path (optional: hafoe will use the default installation path in the user's directory)
\n\t-cdhitest2d\tcd-hit-est-2d path (optional: hafoe will use the default installation path in the user's directory)
\n\t-clustalo\tclustalo path (optional: hafoe will use the default installation path in the user's directory)
\n\t-rlib\tpath to the directory where newly installed R libraries should be stored (optional: hafoe will create rlib directory in the output directory by default)
\n\t-readlength\tfragment size used for neighbor-aware serotype identification (optional: default is 100)
\n\t-stepsize\tthe distance between consecutive fragments (optional: the default is 10)
\n\t-vd_criterion\toption used for filtering out parental serotypes of low quality, values include sum, avg (optional: default is sum)
\n\t-title_of_the_run\ttitle of the run (optional: chimericlib filename is used as default)
\n"



#######################################################################
############### 	argument check		 ######################
#######################################################################

if [[ $# == 0 ]]; then
  echo -e $usage
  exit 1	
fi 



i=1
let n=$#+1
out="hafoe_out"
while [ $i -lt $n ]; do
	declare arg=${!i}
	case $arg in 
		"-parentlib")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <fa> argument not specified after \"-parentlib\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare parentlib=${!j}			
			fi		
		;;	

		"-chimericlib")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <csv> argument not specified after \"-chimericlib\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare chimericlib=${!j}			
			fi
		;;	

		"-enrichedlib1")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <fq1> argument not specified after \"-enrichedlib1\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare enrichedlib1=${!j}			
			fi
		;;	
		"-enrichedlib2")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <fq2> argument not specified after \"-enrichedlib2\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare enrichedlib2=${!j}							
			fi
		;;
		"-exploreout")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <path_to_output_folder> argument not specified after \"-exploreout\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare exploreout=${!j}							
			fi
		;;
		"-o")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <o> argument not specified after \"-o\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare out=${!j}							
			fi
		;;
		"-samtools")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-samtools> argument not specified after \"samtools\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare samtools=${!j}	
			fi
		;;	
		"-bowal")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-bowal> argument not specified after \"bowal\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare bowtie_align=${!j}	
			fi
		;;	
		"-bowb")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-bowb> argument not specified after \"bowb\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare bowtie_build=${!j}	
			fi
		;;
		"-cdhitest")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-cdhitest> argument not specified after \"cdhitest\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare cd_hit_est=${!j}	
			fi
		;;
		"-cdhitest2d")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-cdhitest2d> argument not specified after \"cdhitest2d\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare cd_hit_est_2d=${!j}	
			fi
		;;
		"-clustalo")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-clustalo> argument not specified after \"clustalo\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare clustalo=${!j}	
			fi
		;;
		"-readlength")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-readlength> argument not specified after \"readlength\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare readlength=${!j}	
			fi
		;;
		"-stepsize")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-stepsize> argument not specified after \"stepsize\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare stepsize=${!j}	
			fi
		;;
		"-vd_criterion")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-vd_criterion> argument not specified after \"vd_criterion\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare vd_criterion=${!j}	
			fi
		;;
		"-title_of_the_run")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-title_of_the_run> argument not specified after \"title_of_the_run\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare title_of_the_run=${!j}	
			fi
		;;
		"-rlib")
			if [ $[$n-$i] -eq 1 ]; then
				echo -e "$(tput setaf 1;) \nError: the <-rlib> argument not specified after \"rlib\"$(tput sgr0)"
				echo -e $usage
				exit 1	
			else
				let j=$i+1
				declare rlib=${!j}	
			fi
		;;
		"--explore")
			declare explore=true
		;;
		"--identify")
			declare identify=true
		;;
		#doesn't work
		"-*|--*")
      echo -e "$(tput setaf 1;) \nError: Unknown option $arg $(tput sgr0)"
      echo -e $usage
      exit 1
    ;;
		
	esac

	let i+=1
done


############### end of argument check ######################



#######################################################################
############### 	setup test 		#######################
#######################################################################

if [ -z $bowtie_build ]; then
	bowtie_build=$(which bowtie2-build)
#	bowtie_build="$computel_dir/bin/bowtie2-build"
	if [ -z $bowtie_build ]; then
		echo -e "$(tput setaf 1;) \nerror: bowtie2-build is not installed on your system.\nEither specify it with the -bowb option or make sure it's in the system path $(tput sgr0)"
		exit 1
	else
		bowtie_build=`readlink -m $bowtie_build`
		if [ -z $bowtie_build ]; then
			echo -e "$(tput setaf 1;) \nerror: Could not find the default location of bowtie_build on your system.\nEither specify it with the -bowb option or make sure it's in system path $(tput sgr0)"
			exit 1
		fi
	fi
fi


if [ -z $bowtie_align ]; then
    bowtie_align=$(which bowtie2-align)
    if [ -z $bowtie_align ]; then
	bowtie_align=$(which bowtie2)
	if [ -z $bowtie_align ]; then
		echo -e "$(tput setaf 1;) \nerror: bowtie2-align or bowtie2 is not installed on your system.\nEither specify it with the -bowal option or make sure it's in the system path $(tput sgr0)"
		exit 1
	else
		bowtie_align=`readlink -m $bowtie_align`
		if [ -z $bowtie_align ]; then
			echo -e "$(tput setaf 1;) \nerror: Could not find the default location of bowtie_align on your system.\nEither specify it with the -bowal option or make sure it's in system path $(tput sgr0)"
			exit 1
		fi
	fi
    fi
    #	bowtie_align="$computel_dir/bin/bowtie2-align"
fi


if [ -z $samtools ]; then
	samtools=`which samtools`
	if [ -z $samtools ]; then
		echo -e "$(tput setaf 1;) \nerror: samtools is not installed on your system.\nEither specify it with the -samtools option or make sure it's in system path $(tput sgr0)"
		exit 1
	else
		samtools=`readlink -m $samtools`
		if [ -z $samtools ]; then
			echo -e "$(tput setaf 1;) \nerror: Could not find the default location of samtools on your system.\nEither specify it with the -samtools option or make sure it's in system path $(tput sgr0)"
			exit 1
		fi
	fi
fi

if [ -z $cd_hit_est ]; then
	cd_hit_est=`which cd-hit-est`
	if [ -z $cd_hit_est ]; then
		echo -e "$(tput setaf 1;) \nerror: cd-hit-est is not installed on your system.\nEither specify it with the -cdhitest option or make sure it's in system path $(tput sgr0)"
		exit 1
	else
		cd_hit_est=`readlink -m $cd_hit_est`
		if [ -z $cd_hit_est ]; then
			echo -e "$(tput setaf 1;) \nerror: Could not find the default location of cd-hit-est on your system.\nEither specify it with the -cdhitest option or make sure it's in system path $(tput sgr0)"
			exit 1
		fi
	fi
fi


if [ -z $cd_hit_est_2d ]; then
	cd_hit_est_2d=`which cd-hit-est-2d`
	if [ -z $cd_hit_est_2d ]; then
		echo -e "$(tput setaf 1;) \nerror: cd-hit-est-2d is not installed on your system.\nEither specify it with the -cdhitest2d option or make sure it's in system path $(tput sgr0)"
		exit 1
	else
		cd_hit_est_2d=`readlink -m $cd_hit_est_2d`
		if [ -z $cd_hit_est_2d ]; then
			echo -e "$(tput setaf 1;) \nerror: Could not find the default location of cd-hit-est-2d on your system.\nEither specify it with the -cdhitest2d option or make sure it's in system path $(tput sgr0)"
			exit 1
		fi
	fi
fi


if [ -z $clustalo ]; then
	clustalo=`which clustalo`
	if [ -z $clustalo ]; then
		echo -e "$(tput setaf 1;) \nerror: clustalo is not installed on your system.\nEither specify it with the -clustalo option or make sure it's in system path $(tput sgr0)"
		exit 1
	else
		clustalo=`readlink -m $clustalo`
		if [ -z $clustalo ]; then
			echo -e "$(tput setaf 1;) \nerror: Could not find the default location of clustalo on your system.\nEither specify it with the -clustalo option or make sure it's in system path $(tput sgr0)"
			exit 1
		fi
	fi
fi

############### end of setup test ######################




#######################################################################
############### 	input files check		###############
#######################################################################

if [ -z $explore ]; then
	explore=false
fi
if [ -z $identify ]; then
	identify=false
fi


if [ "$explore" = true ] && [ "$identify" != true ] ; then
  if [ -z $parentlib ]; then
  	echo -e "$(tput setaf 1;) \nError:\tthe <fa> argument is required but is not specified $(tput sgr0)"
  	echo -e $usage
  	exit 1
  fi
  
  if [ -z $chimericlib ]; then
  	echo -e "$(tput setaf 1;) \nError:\tthe <csv> argument is required but is not specified $(tput sgr0)"
  	echo -e $usage
  	exit 1
  fi
  
elif [ "$explore" != true ] && [ "$identify" == true ] ; then
  echo 'Be careful not to fall off!'
  if [ -z $exploreout ]; then
  	echo -e "$(tput setaf 1;) \nError:\tthe <path_to_output_folder> argument is required but is not specified $(tput sgr0)"
  	echo -e $usage
  	exit 1
  fi
  
elif [ "$explore" == true ] && [ "$identify" == true ] ; then
  if [ -z $parentlib ]; then
  	echo -e "$(tput setaf 1;) \nError:\tthe <fa> argument is required but is not specified $(tput sgr0)"
  	echo -e $usage
  	exit 1
  fi
  
  if [ -z $chimericlib ]; then
  	echo -e "$(tput setaf 1;) \nError:\tthe <csv> argument is required but is not specified $(tput sgr0)"
  	echo -e $usage
  	exit 1
  fi
  
  if [ -z $enrichedlib1 ]; then
  	echo -e "$(tput setaf 1;) \nError:\tthe <fq1> argument is required but is not specified $(tput sgr0)"
  	echo -e $usage
  	exit 1
  fi
  
  #one enriched is enough
  # if [ -z $enrichedlib2 ]; then
  # 	echo -e "$(tput setaf 1;) \nError:\tthe <fq2> argument is required but is not specified $(tput sgr0)"
  # 	echo -e $usage
  # 	exit 1
  # fi
  
else
  echo -e "$(tput setaf 1;) \nError:\tat least one of the --explore, --identify options should be specified $(tput sgr0)"
	echo -e $usage
	exit 1
fi

#check validity of files

if [ ! -z $parentlib ]; then
	if [ ! -f $parentlib ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe file <fa> $parentlib does not exist $(tput sgr0)"
		echo -e $usage
		exit 1
	fi
fi

if [ ! -z $chimericlib ]; then
	if [ ! -f $chimericlib ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe file <csv> $chimericlib does not exist $(tput sgr0)"
		echo -e $usage
		exit 1
	fi
fi

if [ ! -z $enrichedlib1 ]; then
	if [ ! -d $enrichedlib1 ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe directory <fq1> $enrichedlib1 is not a valid directory. Please provide a valid directory. $(tput sgr0)"
		echo -e $usage
		exit 1
	fi
fi

if [ ! -z $enrichedlib2 ]; then
	if [ ! -d $enrichedlib2 ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe directory <fq2> $enrichedlib2 is not a valid directory. Please provide a valid directory. $(tput sgr0)"
		echo -e $usage
		exit 1
	fi
fi


if [ ! -z $exploreout ]; then
	if [ ! -d $exploreout ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe directory <path_to_output_folder> $exploreout is not a valid directory. Please provide a valid directory. $(tput sgr0)"
		echo -e $usage
		exit 1
	fi
fi

if [ ! -z $rlib ]; then
	if [ ! -d $rlib ]; then
		echo -e "$(tput setaf 1;) \nError:\tthe directory $rlib is not a valid directory. Please provide a valid directory. $(tput sgr0)"
		echo -e $usage
		exit 1
	fi
fi




# try not to overwrite the output file
if [ -f $out ]; then
	echo -e "$(tput setaf 1;)\nError:$(tput sgr0;)\tthe output directory $out is a regular file. Please provide a valid directory."
	exit 1
fi


if [ -d $out ]; then
	#echo -e "$(tput setaf 3;)\nWarning:$(tput sgr0;)\tthe output directory $out already exists. Do you want to replace it? (y/n) "
	echo -e "$(tput setaf 3;)\nWarning:$(tput sgr0;)\tthe output directory $out already exists. Files will be replaced. "
	rm -r $out
	mkdir $out
	#ans=""
	#read ans
	ans="y"
	while [ ! "$ans" == "y" ]; do
		if [ "$ans" == "n" ]; then
			echo "hafoe cancelled"
			exit 1
		fi
		read ans
	done
else
	mkdir $out
fi

	############### end of input files check ###################	
#######################################################################	
############	 setting configuration parameters	###############	
#######################################################################	
program_dir=`dirname $0`	
if [ ! -z $enrichedlib1 ]; then	
	if [ ! -z $enrichedlib2 ]; then	
		enriched="$enrichedlib1,$enrichedlib2"	
	else enriched="$enrichedlib1"	
	fi	
else	
	enriched="none"	
fi	
if [ -z $exploreout ]; then	
	exploreout="none"	
fi	
if [ -z $parentlib ]; then	
	parentlib="none"	
fi	
if [ -z $chimericlib ]; then	
	chimericlib="none"	
fi
if [ -z $readlength ]; then	
	readlength=100	
fi
if [ -z $stepsize ]; then	
	stepsize=10
fi
if [ -z $vd_criterion ]; then	
	vd_criterion="sum"
fi
if [ -z $title_of_the_run ]; then	
  title_of_the_run_=$(basename $chimericlib)
	title_of_the_run="${title_of_the_run_%.*}"
fi
if [ -z $rlib ]; then	
  rlib=$out/rlib
  mkdir $out/rlib
fi

############### end of setting configuration parameters ######################



#######################################################################
##########	generating config file and running	###############
#######################################################################

if [ ! -d $out/config ]; then
  mkdir $out/config
fi


# generate config file
cat > "$out/config/config_file" << bzez
scripts.dir $program_dir/src
bowtie.build.path $bowtie_build
bowtie.align.path $bowtie_align
samtools.path $samtools
cd_hit_est.path $cd_hit_est
cd_hit_est_2d.path $cd_hit_est_2d
clustalo.path $clustalo
explore $explore
identify $identify
parent $parentlib
chimera $chimericlib
enriched $enriched
explore.out $exploreout
output.dir $out
read_length $readlength
step_size $stepsize
vd_criterion $vd_criterion
title_of_the_run $title_of_the_run
rlib $rlib
bzez


chmod +x $program_dir/src/*.sh
Rscript $program_dir/src/validate.options.R $out/config/config_file

#After running the validate.options.R and pipline.R, 
#run viz/generate_plots.py to generate bokeh reports
python3 $program_dir/viz/generate_plots.py $out
