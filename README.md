# _hafoe_ 
_hafoe_ is a command-line-based tool for the automated exploratory analysis of AAV chimeric libraries and identification of enriched variants in desired tissues using long-read sequencing datasets.

### Operating Systems
_hafoe_ works with Unix operating system (tested for Ubuntu Linux).

### Installation
Download and uncompress the _hafoe_ package in a local directory. Then make hafoe.sh file executable by running `chmod +x hafoe.sh`.

### Usage
The common _hafoe_ usage is:

```
./hafoe.sh \
    --explore \
    --identify \
    -parentlib <path_to_parental_sequences_fasta_file> \
    -chimericlib <path_to_chimeric_sequences_csv/fastq_file> \
    -enrichedlib1 <path_to_dir_with_enriched_sequences_fastq_file(s)> \
    -o <output_directory> example/hafoe_out \
    -title_of_the_run <title> 
```

### Required arguments
_hafoe_ works with one of the --explore, --identify options or both

--explore option should be specified for exploratory analysis of chimeric sequences. When using this option, the required arguments are: 
```
-parentlib the full path to the fasta file of parental sequences
-chimericlib the full path to the csv file of chimeric sequences and their abundances or fastq file of chimeric sequences
```

--identify option should be specified for identification of novel tissue-specific variants. When using this option, the required argument is: 
```
-enrichedlib1 the full path to directory containing fastq file(s) of sequences obtained after enrichment in one or more tissue samples
```

When using --identify option without --explore, the required argument is:
```
-exploreout the full path to output directory generated by running _hafoe_ with only --explore option
```

### Additional arguments
Run `./hafoe.sh` to see additional arguments.

### Test run
To test how _hafoe_ works navigate to hafoe/ directory and run `example/run_hafoe.sh`. All required parameters are specified, however, you will need to change -cdhitest, -cdhitest2d, -rlib options to point to corresponding tools in your system.

A successful test run should generate the example/hafoe_out/ output directory with output files and reports.

### Output
_hafoe_'s main outputs are files/, reports/ directories in output directory.

### System requirements
The following software and packages should be pre-installed in your system: 

R (v4.1.3) with the following packages:
dplyr (v1.0.9),
ORFik (v1.12.13),
plotly (v4.10.0),
ggplot2 (v3.3.6),
gplots (v3.1.3),
microseq (v2.1.4),
Biostrings (v2.60.2),
string (v1.4.0),
cowplot (v1.1.1),
seqinr (v4.2.8),

Python (v3.9.5) with the following packages: 
numpy (v1.24.3),
pandas (v1.3.4),
Bio (v1.5.3),
bokeh (v2.4.3),
seaborn (v0.12.2),
selenium (v4.10.0),

Bowtie 2 (v2.4.2), 

CD-HIT (v4.8.1), 

Clustal Omega (v1.2.4), 

SAMtools (v1.9).
 
