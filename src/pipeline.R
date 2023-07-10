
rdata_name <- paste0("./", output.dir, ".RData")

# Appending the log file
#######################################################################################
#######################################################################################

#log_file <- file(file.path(output.dir, "log/hafoe.log"), open = "wt")
sink(log_file, append = T, type = "output") 
sink(log_file, append = T, type = "message")

# Checking and installing dependencies
#######################################################################################
#######################################################################################

cat("\nChecking and installing dependencies\n")
cat("===================================================================\n\n")


# If a package is installed, it will be loaded. If any
# are not, the missing package(s) will be installed
# and then loaded.

# First specify the packages of interest
packages = c("dplyr", "ORFik", "plotly", "ggplot2", "gplots",
             "microseq", "Biostrings", "stringr", "cowplot", "seqinr",
             "reshape2", "ggrepel", "fmsb", "RColorBrewer")

.libPaths(rlib)

# Now load or install&load all
if (!suppressPackageStartupMessages(require("BiocManager", quietly = TRUE))){
  install.packages("BiocManager", lib = rlib)
}

if (!suppressPackageStartupMessages(require("ORFik", character.only = TRUE, quietly = T))) {
  BiocManager::install("ORFik",  lib = rlib)
  suppressPackageStartupMessages(library("ORFik", character.only = TRUE, quietly = T))
}

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = T))) {
      install.packages(x, dependencies = TRUE, lib = rlib)
      suppressPackageStartupMessages(library(x, character.only = TRUE, quietly = T))
    }
  }
)


cat(paste0("Following R packages were installed in ", rlib, " directory and/or called:\n"))
cat(paste0("\tdplyr: version ", packageVersion("dplyr"), "\n"))
cat(paste0("\tORFik: version ", packageVersion("ORFik"), "\n"))
cat(paste0("\tplotly: version ", packageVersion("plotly"), "\n"))
cat(paste0("\tggplot2: version ", packageVersion("ggplot2"), "\n"))
cat(paste0("\tgplots: version ", packageVersion("gplots"), "\n"))
cat(paste0("\tmicroseq: version ", packageVersion("microseq"), "\n"))
cat(paste0("\tBiostrings: version ", packageVersion("Biostrings"), "\n"))
cat(paste0("\tstringr: version ", packageVersion("stringr"), "\n"))
cat(paste0("\tcowplot: version ", packageVersion("cowplot"), "\n"))
cat(paste0("\tseqinr: version ", packageVersion("seqinr"), "\n"))

source(file.path(scripts.dir, "functions.R"))


# Check what part of the program to run
#######################################################################################
#######################################################################################

if (explore) {


  # Read, pre-process the input files
  #######################################################################################
  #######################################################################################

  parent_serotypes <- microseq::readFasta(parent.path)

  if (grepl(".csv", chimeric.lib.path)){
    chimeric_library <- read.csv(chimeric.lib.path) #should contain columns: X (name), Sequence, Count      
    if (any(chimeric_library$X == "")){
      # Handle missing names
      index <- c(1:nrow(chimeric_library))
      chimeric_library[, "index"] <- index
      chimeric_library[, "X"] = paste0("AAV.", 100000 + as.numeric(chimeric_library$index))
    }
    
    chimeric_library_fa <- chimeric_library[c('X', 'Sequence')]
    colnames(chimeric_library_fa) <- c("Header", "Sequence")
    chimeric_library_fa.path <- paste0(sub('\\..[^\\.]*$', '', chimeric.lib.path), ".fasta")
    microseq::writeFasta(chimeric_library_fa, chimeric_library_fa.path)
  } else if (grepl(".fastq", chimeric.lib.path)){
    #deduplicate fastq sequneces to get chimeric_library dataframe
    chimeric_library_fq <- microseq::readFastq(chimeric.lib.path)    
    chimeric_library <- deduplicate(file = chimeric_library_fq, 
                                    out_path = paste0(output.dir, "/files/Chimeric_lib_unique_seq.csv")) 
    
    chimeric_library_fa <- chimeric_library[c('X', 'Sequence')]
    colnames(chimeric_library_fa) <- c("Header", "Sequence")
    chimeric_library_fa.path <- paste0(sub('\\..[^\\.]*$', '', sub('\\..[^\\.]*$', '', chimeric.lib.path)), ".fasta")
    microseq::writeFasta(chimeric_library_fa, chimeric_library_fa.path)
  }


  # Variant length distribution
  #######################################################################################

  # p1 <- plot.histogram(chimeric_library_fa$Sequence,
  #                      "Chimeric",
  #                      "variant",
  #                      bin_size = 100)
  # p2 <- plot.histogram(enriched1$Sequence,
  #                      "Enriched1",
  #                      "variant",
  #                      bin_size = 100)
  # p3 <- plot.histogram(enriched2$Sequence,
  #                      "Enriched2",
  #                      "variant",
  #                      bin_size = 100)


  # Filter by ORF length >= min_orf (600 aa) and sequence length <= max_seq(3000 nt)
  # ORF, Variant length distribution after filtering
  #######################################################################################
  #######################################################################################
  dir.create(file.path(output.dir, "files"), showWarnings = F)
  ### save output paths in config/output_paths
  output_paths <- data.frame(file_name = NA,
                             path = NA)

  cat("\n\nFiltering the chimeric library by ORF and original sequence lengths\n")
  cat("===================================================================\n\n")

  save.image(rdata_name)  

  # Chimeric library
  out_chim <- filter.by.orf(file = chimeric_library_fa, 
                            file_path = chimeric_library_fa.path,
                            library_name = "Chimeric", 
                            output_paths = output_paths)
  fasta_orf_filtered_chim <- out_chim[[1]]
  fasta_orf_seq_filtered_chim <- out_chim[[2]]

  output_paths <- rbind(output_paths, 
                        c('final_chimeric_orf_fasta', file.path(output.dir, "files", paste0("Chimeric", "_ORF.fasta"))))
  output_paths <- rbind(output_paths,
                        c("chimeric_orf_summary", file.path(output.dir, "files", "Chimeric_orf_summary.csv")))
  
  p4 <- plot.histogram(fasta_orf_filtered_chim$Sequence,
                       "Chimeric",
                       "variant",
                       additional_info = "\nafter filtering")
  p5 <- plot.histogram(fasta_orf_seq_filtered_chim$Sequence,
                       "Chimeric",
                       "ORF",
                       additional_info = "\nafter filtering")


  # Clustering filtered chimeric library
  #######################################################################################
  #######################################################################################

  cat("\n\nClustering the filtered chimeric library\n")
  cat("===================================================================\n\n")
  cat("Program: cd-hit-est\n")


  dir.create(file.path(output.dir, "files/clstr_chimeric_lib"), showWarnings = F)

  o <- system(paste0(file.path(scripts.dir, "clustering.sh"),
                " -h ", cd_hit_est.path,
                " -i ", file.path(output.dir, "files", paste0("Chimeric", "_ORF.fasta")),
                " -o ", file.path(output.dir, "files/clstr_chimeric_lib"),
                " -p ", "clstr_chim.fasta ",
                " -l ", file.path(output.dir, "log/cd_hit_chimeric.log"),
                " -c ", 0.90,  
                " -n ", 9,
                " -g ", 0,
                " -a ", 0.90),
          intern = TRUE) #to capture the output of the command as an R character vector
  cat(o, sep = "\n")


  cluster_sizes <- read.csv(file.path(output.dir, "files/clstr_chimeric_lib/cluster_sizes.csv"), header = F)


  cat(paste0("\nNumber of clusters: ", nrow(cluster_sizes), "\n"))
  cat(paste0("Cluster size summary statistics\n"))
  print(summary(cluster_sizes$V1))



  # Choose representatives
  #######################################################################################
  #######################################################################################

  cat("\n\nChoosing new representative sequences for the clusters\n")
  cat("===================================================================\n\n")

  cluster_members <- choose.new.representatives(sizes_path = file.path(output.dir, "files/clstr_chimeric_lib/cluster_sizes.csv"),
                                               members_path = file.path(output.dir, "files/clstr_chimeric_lib/members_ordered.csv"),
                                               chimeric_library = chimeric_library,
                                               fasta_orf_seq_filtered_chim = fasta_orf_seq_filtered_chim,
                                               fasta_path = file.path(output.dir, "files/clstr_chimeric_lib/clstr_chim_new.fasta"),
                                               csv_path = file.path(output.dir, "files/clstr_chimeric_lib/representatives.csv"))
  
  output_paths <- rbind(output_paths,
                        c("chimeric_lib_clustering_info", file.path(output.dir, "files/clstr_chimeric_lib/clustering_info.csv")))
  
  
  cat("\n\nPlotting the cluster size distribution\n")
  cat("===================================================================\n\n")

  dir.create(file.path(output.dir, "reports"), showWarnings = F)

  cat(paste0("Output path: ", file.path(output.dir, "reports/cluster_size_distr_chimeric_lib.pdf"), "\n"))

  graphics.off()
  pdf(file.path(output.dir, "reports/cluster_size_distr_chimeric_lib.pdf"), width=8, height=5)
  plot.cluster.size.distribution(cluster_members = cluster_members,
                                 size_thresh = 45,
                                 similarity_thresh = 90,
                                 library_name = "chimeric library")
  dev.off()


  #Index parental AAV file
  #######################################################################################
  #######################################################################################

  cat("\n\nBuilding index on parental library\n")
  cat("===================================================================\n\n")


  dir.create(file.path(output.dir, "files/parent_index"), showWarnings = F)

  o <- system(paste0(file.path(scripts.dir, "index.sh"),
                      " -b ", bowtie.build.path,
                      " -i ", parent.path,
                      " -o ", file.path(output.dir, "files/parent_index"),
                      " -p ", "AAV",
                      " -l ", file.path(output.dir, "log/bowtie_build.log")),
              intern = TRUE)
  cat(o, sep = "\n")



  # Variant description of chimeric library representatives
  #######################################################################################
  #######################################################################################

  cat("\n\nDoing variant description on chimeric library representatives\n")
  cat("===================================================================\n")

  cat("\n\nMaking fastq files for chimeric library representatives\n")
  cat("-------------------------------------------------------------------\n\n")

  #make fastq files for query sequences
  make.fastq.files(read_length = read_length,
                   fasta.path = file.path(output.dir, "files/clstr_chimeric_lib/clstr_chim_new.fasta"),
                   library_name = "chimeric_lib_representatives",
                   step_size = step_size)
  
  cat("\n\nAligning chimeric library representatives on parental index\n")
  cat("-------------------------------------------------------------------\n\n")

  #do alignment
  o <- system(paste0(file.path(scripts.dir, "align_strict.sh"),
                      " -b ", bowtie.align.path,
                      " -d ", file.path(output.dir, "files/variant_description/chimeric_lib_representatives"),
                      " -r ", file.path(output.dir, "files/parent_index", "AAV"),
                      " -s ", samtools.path,
                      " -i ", "fastq",
                      " -l ", file.path(output.dir, "log/bowtie_align_representatives.log"),
                      " -L ", 30,
                      " -D ", 2 ),
              intern = T)
  cat(o, sep = "\n")


  cat("\n\nGetting required files for variant description\n")
  cat("-------------------------------------------------------------------\n\n")

  #get description files
  o <- system(paste0(file.path(scripts.dir, "get_description_files.sh"),
                     " -d ", file.path(output.dir, "files/variant_description/chimeric_lib_representatives"),
                     " -s ", samtools.path),
              intern = T)
  cat(o, sep = "\n")


  #Variant description


  #list mapping AAV serotypes' accession codes to number
  serotype_name <- c(
    "NC_002077.1",
    "NC_001401.2",  #old "ENA|J01901|J01901.1"
    "ENA|U48704|U48704.1" ,
    "NC_001829.1",
    "NC_006152.1",
    "ENA|AF028704|AF028704.1",
    "ENA|AF513851|AF513851.1",
    "ENA|AAN03857|AAN03857.1",
    "ENA|AAS99264|AAS99264.1",
    "ENA|AAT46337|AAT46337.1",
    "AY631966.1",
    "DQ813647.1",
    "EU285562.1",
    "AY242997.1", #AAVrh8
    "AY243015.1", #AAVrh10
    "AY243003.1") #AAVrh32
  serotype_num <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
  serotypes_df <- data.frame(serotype_name, serotype_num, stringsAsFactors = FALSE)


  
  col = c("#D3D3D3", "#A6CEE3", "#1F78B4", "#B2DF8A", "#555555", "#33A02C",
          "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
          "#FFFF99", "#B15928", "yellow", "#999999", "#a70000", "black", "white")
  legend_name <- c("no alignment", "AAV1", "AAV2", "AAV3", "AAV4", "AAV5",
                   "AAV6", "AAV7", "AAV8", "AAV9", "AAV10", "AAV11", "AAV12", 
                   "AAV13", "AAVrh8", "AAVrh10", "AAVrh32", "multiple alignment", "gap")
  col_df <- data.frame(col = col)
  rownames(col_df) <- legend_name
  
  
  
  save.image(rdata_name)

  cat("\n\nDoing neighbor-aware serotype identification\n")
  cat("-------------------------------------------------------------------\n\n")

  neighbor.joining(parents_df = serotypes_df,
                   path_fasta = file.path(output.dir, "files/clstr_chimeric_lib/clstr_chim_new.fasta"),
                   path_mapped = file.path(output.dir, "files/variant_description/chimeric_lib_representatives/csv/mapped"),
                   path_unmapped = file.path(output.dir, "files/variant_description/chimeric_lib_representatives/csv/unmapped"),
                   read_length = read_length,
                   library_name = "chimeric_lib_representatives",
                   step_size = step_size, 
                   vd_criterion = vd_criterion)

  output_paths <- rbind(output_paths,
                        c("representatives_variant_description", 
                          file.path(output.dir, 
                                    "files/variant_description/chimeric_lib_representatives/chimeric_lib_representatives_variant_description.csv")))
  

  #Plot variant description

  reps_nj_mat <- read.table(file.path(output.dir,
                                      "files/variant_description/chimeric_lib_representatives/chimeric_lib_representatives_variant_description.csv"))
  rownames(reps_nj_mat) <- reps_nj_mat[,1]
  reps_nj_mat[,1] <- NULL
  reps_nj_mat <- as.matrix(reps_nj_mat)
  colnames(reps_nj_mat) <- NULL
  

  serotypes_freq <- get.frequency.table(chimeric_library,
                                        nj_matrix = reps_nj_mat,
                                        cluster_members,
                                        output_path = file.path(output.dir, "files"))
  
  output_paths <- rbind(output_paths,
                        c("chimeric_lib_serotype_frequency", file.path(output.dir, "files", "chimeric_lib_serotype_frequency.csv")))
         
  output_paths <- rbind(output_paths,
                        c("chimeric_lib_representatives_counts", file.path(output.dir, "files", "chimeric_lib_rep_counts.csv")))
  
  cat("\n\nPlotting the cluster abundance\n")
  cat("===================================================================\n\n")

  cat(paste0("Output path: ", file.path(output.dir, "reports/cluster_abundance_chimeric_lib.pdf"), "\n"))

  graphics.off()
  pdf(file.path(output.dir, "reports/cluster_abundance_chimeric_lib.pdf"), width=8, height=5)
  plot.cluster.abundance(file_path = file.path(output.dir, "files/chimeric_lib_rep_counts.csv"),
                         size_thresh = 45,
                         similarity_thresh = 90,
                         library_name = "chimeric library")
  dev.off()

  cat("\n\nPlotting variant description of chimeric library representatives\n")
  cat("===================================================================\n\n")
  cat(paste0("Output path: ", file.path(output.dir, "reports/variant_description_chimeric_rep.pdf"), "\n"))

  graphics.off()
  pdf(file.path(output.dir, "reports/variant_description_chimeric_rep.pdf"), width=8, height=7)
  plot.variant.description(matrix = reps_nj_mat, 
                           col_df = col_df,
                           library_name = "chimeric library representatives\n")
  dev.off()
  
  cat("\n\nMultiple sequence alignment of top 20 chimeric library representatives\n")
  cat("===================================================================\n\n")
  cat(paste0("Output path: ", file.path(output.dir, "files/top20_rep_msa.clustal_num"), "\n"))

  cluster_chim_new <- microseq::readFasta(file.path(output.dir, "files/clstr_chimeric_lib/clstr_chim_new.fasta"))
  counts_chim <- read.csv(file.path(output.dir, "files/chimeric_lib_rep_counts.csv"))
  
  top20_rep_fa <- cluster_chim_new[gsub("\\/.*","",cluster_chim_new$Header) %in% counts_chim[1:20,]$Representative,]
  top20_rep_fa$Header <- gsub("\\/.*","",top20_rep_fa$Header)
  microseq::writeFasta(top20_rep_fa, file.path(output.dir, "files/top20_representatives.fasta"))
  
  o <- system(paste0(file.path(scripts.dir, "msa.sh"),
                     " -c ", clustalo.path,
                     " -i ", file.path(output.dir, "files/top20_representatives.fasta"),
                     " -o ", file.path(output.dir, "files/top20_rep_msa.clustal_num")),
              intern = TRUE) 
  cat(o, sep = "\n\n")
  
  cat("\n\nGet conserved and gap regions from multiple sequence alignment\n")
  cat("===================================================================\n\n")
  cat(paste0("Output path: ", file.path(output.dir, "files/"), "\n"))

  out <- get.conserved.positions(aln_file_path = file.path(output.dir, "files/top20_rep_msa.clustal_num"))
  identity_ranges_top20 <- out[[1]]
  alignment_top20 <- out[[2]]

  
  matrix_nt_top20 <- add.gap.info(alignment = alignment_top20, 
                                  nj_matrix = reps_nj_mat, 
                                  step_size = step_size, 
                                  nj_by_nt = T) 

  write.table(matrix_nt_top20, 
              file.path(output.dir, 
                        "files/variant_description", 
                        "chimeric_lib_representatives", 
                        paste0("chimeric_lib_representatives", "_variant_description_top20_msa.csv")), 
              #sep = ",",
              quote = F, 
              col.names = F)
  
  output_paths <- rbind(output_paths,
                        c("representatives_variant_description_top20", 
                          file.path(output.dir, 
                                    "files/variant_description/chimeric_lib_representatives/chimeric_lib_representatives_variant_description_top20_msa.csv")))
  
  
  
  cat("\n\nPlot variant description of top20 representatives\nwith multiple sequence alignment gaps\n")
  cat("===================================================================\n\n")
  cat(paste0("Output path: ", file.path(output.dir, "reports/variant_description_top20_msa.pdf"), "\n"))

  graphics.off()
  pdf(file.path(output.dir, "reports/variant_description_top20_msa.pdf"), width=8, height=6)
  plot.variant.description(matrix = matrix_nt_top20, 
                           col_df = col_df,
                           library_name = "chimeric library top 20 representatives\n")
  dev.off()
  
  
  
  
  #Plot VP1,2,3 positions
  # graphics.off()
  # pdf(file.path(output.dir, "reports/viral_proteins_positions.pdf"), width=8, height=1.5)
  # plot.vp.positions()
  # dev.off()
  
  
  cat("\n\nPlotting distribution of AAV serotypes in chimeric library\nbased on variant description of representatives and their\nabundance the chimeric library\n")
  cat("===================================================================\n\n")

  cat(paste0("Output path: ", file.path(output.dir, "reports/serotype_distribution_chimeric_lib.pdf"), "\n"))

  graphics.off()
  pdf(file.path(output.dir, "reports/serotype_distribution_chimeric_lib.pdf"), width=8, height=5)
  plot.serotype.frequency(serotypes_freq = serotypes_freq,
                          col_df = col_df,
                          library_name = "chimeric library")
  dev.off()


  get.reps.nj.matrix.nt(output_path = file.path(output.dir, "files"),
                        path_fasta = file.path(output.dir, "files/clstr_chimeric_lib/clstr_chim_new.fasta"),
                        nj_matrix = reps_nj_mat,
                        nj_by_nt = T) #with updated neighbor_joining function T

  if(identify) {
    # Use enriched libs to find abundant variants

    # Read, pre-process Enriched1 library files
    #########################################
    #########################################
    
    for (en1.path in en1.paths){
      if (!(is.na(en1.path) || en1.path == '' || is.null(en1.path))) {
        en1_name <- basename(en1.path)
        en1_name <- gsub("\\..*", "", en1_name) 
        
        dir.create(file.path(output.dir, "files", paste0("enriched1_", en1_name)), showWarnings = F)
        
        cat(paste0("\n\nFiltering the enrichment cycle 1 library ", en1_name,"\nby ORF and original sequence lengths\n"))
        cat("===================================================================\n\n")
        
        enriched1 <- microseq::readFastq(en1.path)
        #enriched1$Header <- gsub("/", ".", enriched1$Header)
        out_e1 <- filter.by.orf(file = enriched1, 
                                file_path = en1.path,
                                library_name = paste0("Enriched1_", en1_name),
                                format = "fastq",
                                output_paths = output_paths)
        fasta_orf_filtered_e1 <- out_e1[[1]]
        fasta_orf_seq_filtered_e1 <- out_e1[[2]]
        output_paths <- out_e1[[3]]
        
        # output_paths <- rbind(output_paths, 
        #                       c('final_enriched_orf_fasta', file.path(output.dir, "files", paste0("Enriched1", "_ORF.fasta"))))
        # output_paths <- rbind(output_paths,
        #                       c("enriched1_orf_summary", file.path(output.dir, "files", "Enriched1_orf_summary.csv")))
        
        
        p6 <- plot.histogram(fasta_orf_filtered_e1$Sequence,
                             "Enriched1",
                             "variant",
                             additional_info = "\nafter filtering")
        p7 <- plot.histogram(fasta_orf_seq_filtered_e1$Sequence,
                             "Enriched1",
                             "ORF",
                             additional_info = "\nafter filtering")
      }
      
      
      # #Enriched2 library
      # if (!(is.na(en2.path) || en2.path == '' || is.null(en2.path))) {
      #   cat("Filtering the enrichment cycle 2 library\n")
      #   cat("-------------------------------------------------------------------\n")
      #   
      #   enriched2 <- microseq::readFastq(en2.path)
      #   out_e2 <- filter.by.orf(enriched2,
      #                           en2.path,
      #                           library_name = "Enriched2",
      #                           format = "fastq")
      #   fasta_orf_filtered_e2 <- out_e2[[1]]
      #   fasta_orf_seq_filtered_e2 <- out_e2[[2]]
      #   
      #   p8 <- plot.histogram(fasta_orf_filtered_e2$Sequence,
      #                        "Enriched2",
      #                        "variant",
      #                        additional_info = "\nafter filtering")
      #   p9 <- plot.histogram(fasta_orf_seq_filtered_e2$Sequence,
      #                        "Enriched2",
      #                        "ORF",
      #                        additional_info = "\nafter filtering")
      #   
      # }
      
      # Clustering enriched library sequences with chimeric library representative sequences
      #######################################################################################
      #######################################################################################
  
      cat("\n\nClustering enriched 1 library with chimeric library representatives\n")
      cat("===================================================================\n\n")
      cat("Program: cd-hit-est-2d\n")
  
      dir.create(file.path(output.dir, "files", paste0("enriched1_", en1_name), "clstr"), showWarnings = F)
  
      o <- system(paste0(file.path(scripts.dir, "clustering_identify.sh"),
                         " -h ", cd_hit_est_2d.path,
                         " -i ", file.path(output.dir, "files/clstr_chimeric_lib/clstr_chim_new.fasta"),
                         " -e ", file.path(output.dir, "files", paste0("Enriched1_", en1_name, "_ORF.fasta")),
                         " -o ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "clstr"),
                         " -p ", "clstr_enriched1.fasta ",
                         " -l ", file.path(output.dir, paste0("log/cd_hit_enriched1_", en1_name, ".log")),
                         " -c ", 0.95,
                         " -n ", 10,
                         " -g ", 0,
                         " -a ", 0.95,
                         " -m ", 8000),
                  intern = TRUE) #to capture the output of the command as an R character vector
      cat(o, sep = "\n")
  
  
      cluster_sizes <- read.csv(file.path(output.dir, "files", paste0("enriched1_", en1_name), "clstr", "cluster_sizes.csv"), header = F)
  
      cat(paste0("\nNumber of clusters: ", nrow(cluster_sizes), "\n"))
      cat(paste0("Cluster size summary statistics\n"))
      print(summary(cluster_sizes$V1))
      
      
      cat("\n\nCreating counts.csv summary table for representative enrichment\n")
      cat("===================================================================\n\n")
      cat(paste0("Output path: ", file.path(output.dir, "files/counts.csv"), "\n"))
      
      get.counts.table(chim_file_path = file.path(output.dir, "files/chimeric_lib_rep_counts.csv"),
                       clstr_enriched1_dir = file.path(output.dir, "files", paste0("enriched1_", en1_name), "clstr"),
                       output_path = file.path(output.dir, "files"), 
                       en1_name = en1_name)
      
      cat("\n\nPlotting top representatives' abundance in enriched 1 library\n")
      cat("===================================================================\n\n")
      
      cat(paste0("Output path: ", file.path(output.dir, "reports", paste0("enriched1_", en1_name), "top_rep_en1_barplot.pdf"), "\n"))
      
      dir.create(file.path(output.dir, "reports", paste0("enriched1_", en1_name)), showWarnings = F)
      
      graphics.off()
      pdf(file.path(output.dir, "reports", paste0("enriched1_", en1_name), "top_rep_en1_barplot.pdf"), width=8, height=5)
      plot.top.reps.in.enrichedlib.ggplot(counts_file_path = file.path(output.dir, "files/counts.csv"),
                                          topn_thresh = 30,
                                          en1_name = en1_name,
                                          library_name = paste0("Enriched1_", en1_name))
      dev.off()
      
      cat("\n\nFilter enriched, reduced representatives based on En1/Chim ratio\n")
      cat("===================================================================\n\n")
      
      cat(paste0("Output path for enriched: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "enriched_representatives.csv"), "\n"))
      cat(paste0("Output path for reduced: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "reduced_representatives.csv"), "\n"))
      cat(paste0("Output path for enriched: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "enriched_representatives.fasta"), "\n"))
      cat(paste0("Output path for reduced: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "reduced_representatives.fasta"), "\n\n"))
      
      save.image(rdata_name)
      
      out <- filter.enriched.reduced.reps(counts_file_path = file.path(output.dir, "files/counts.csv"),
                                          output_path = file.path(output.dir, "files", paste0("enriched1_", en1_name)),
                                          en1_name = en1_name,
                                          reduced_chim_thresh = 0.1,
                                          reduced_ratio_thresh = 1,
                                          enriched_en1_thresh = 0.1,
                                          enriched_ratio_thresh = 1)
      enriched = out[[1]]
      reduced = out[[2]]
      
      
      cat(paste0(nrow(enriched), " enriched sequences found.\n"))
      cat(paste0(nrow(reduced), " reduced sequences found.\n"))
      
      
      cat("\n\nPlot enriched, reduced representatives based on En1/Chim ratio\n")
      cat("===================================================================\n\n")
      cat(paste0("Output path for enriched: ", file.path(output.dir, "reports", paste0("enriched1_", en1_name), "enriched_reps_en1.pdf"), "\n"))
      cat(paste0("Output path for reduced: ", file.path(output.dir, "reports", paste0("enriched1_", en1_name), "reduced_reps_en1.pdf"), "\n"))
      
      graphics.off()
      pdf(file.path(output.dir, "reports", paste0("enriched1_", en1_name), "enriched_reps_en1.pdf"), width=12, height=5)
      plot.enrichment.tiles(enriched = enriched[1:min(nrow(enriched), 20),], 
                            reduced = reduced,
                            type = "Enriched", 
                            en1_name = en1_name) 
                            #axis_text_size = 8)
      dev.off()
      
      reduced <- reduced[order(reduced$Chimeric.Count, decreasing = T),]
      graphics.off()
      pdf(file.path(output.dir, "reports", paste0("enriched1_", en1_name), "reduced_reps_en1.pdf"), width=12, height=5)
      plot.enrichment.tiles(enriched = enriched, 
                            reduced = reduced[1:min(nrow(reduced), 20),],
                            type = "Reduced",
                            en1_name = en1_name)
      dev.off()
      
      cat("\n\nPlotting variant description of enriched, reduced representatives\n")
      cat("===================================================================\n\n")
      cat(paste0("Output path: ", file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_enriched.pdf"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_reduced.pdf"), "\n"))
      
      graphics.off()
      pdf(file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_enriched.pdf"), width=8, height=5)
      plot.variant.description(matrix = reps_nj_mat[enriched[1:min(nrow(enriched), 20),]$Representative,],
                               col_df = col_df,
                               library_name = "enriched representative variants\n")
      dev.off()
      
  
      graphics.off()
      pdf(file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_reduced.pdf"), width=8, height=5)
      plot.variant.description(matrix = reps_nj_mat[reduced[1:min(nrow(reduced), 20),]$Representative,],
                               col_df = col_df,
                               library_name = "reduced representative variants\n")
      dev.off()
  
      cat("\n\nMultiple sequence alignment of enriched, reduced representatives\n")
      cat("===================================================================\n\n")
      cat(paste0("Output path: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "enriched_msa.clustal_num"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "reduced_msa.clustal_num"), "\n\n"))
      
      o <- system(paste0(file.path(scripts.dir, "msa.sh"),
                         " -c ", clustalo.path,
                         " -i ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "enriched_representatives.fasta"),
                         " -o ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "enriched_msa.clustal_num")),
                  intern = TRUE) 
      cat(o, sep = "\n\n")
      
      
      o <- system(paste0(file.path(scripts.dir, "msa.sh"),
                         " -c ", clustalo.path,
                         " -i ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "reduced_representatives.fasta"),
                         " -o ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "reduced_msa.clustal_num")),
                  intern = TRUE) 
      cat(o, sep = "\n")
      
      
      cat("\n\nGet conserved and gap regions from multiple sequence alignment\n")
      cat("===================================================================\n\n")
      cat(paste0("Output path: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "/"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "files", paste0("enriched1_", en1_name), "/"), "\n"))
      
      out <- get.conserved.positions(aln_file_path = file.path(output.dir, "files", paste0("enriched1_", en1_name), "enriched_msa.clustal_num"))
      identity_ranges_enriched <- out[[1]]
      alignment_enriched <- out[[2]]
      
      out <- get.conserved.positions(aln_file_path = file.path(output.dir, "files", paste0("enriched1_", en1_name), "reduced_msa.clustal_num"))
      identity_ranges_reduced <- out[[1]]
      alignment_reduced <- out[[2]]
      
      matrix_nt_enriched <- add.gap.info(alignment = alignment_enriched, 
                                         nj_matrix = reps_nj_mat, 
                                         step_size = step_size, 
                                         nj_by_nt = T) 
      matrix_nt_reduced <- add.gap.info(alignment = alignment_reduced, 
                                        nj_matrix = reps_nj_mat, 
                                        step_size = step_size, 
                                        nj_by_nt = T) 
      
      
      #Use AAV2  Variable regions for plot
      aav2_vr_ranges <- data.frame(VR = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"), 
                              start_aa = c(263, 325, 381, 450, 490, 527, 545, 585, 704), 
                              end_aa = c(265, 330, 384, 466, 503, 532, 556, 596, 713))
      aav2_vr_ranges['start_nt'] <- aav2_vr_ranges$start_aa * 3 - 2    #ATGAGGAGC
      aav2_vr_ranges['end_nt'] <- aav2_vr_ranges$end_aa * 3 
      
      cat("\n\nPlot variant description of enriched and reduced representatives\nwith conserved regions and gaps\n")
      cat("===================================================================\n\n")
      cat(paste0("Output path: ", file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_enriched_conserved.pdf"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_reduced_conserved.pdf"), "\n"))
      
      graphics.off()
      pdf(file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_enriched_conserved.pdf"), width=8, height=5)
      # doesn't recognize identity_ranges
      # plot.variant.description.conserved(matrix = matrix_nt_enriched, identity_ranges = identity_ranges_enriched,
      #                                    col_df = col_df,
      #                                    library_name = "enriched representative variants")
      ####
      matrix = matrix_nt_enriched
      identity_ranges = identity_ranges_enriched
      col_df = col_df
      library_name = "enriched representative variants\n"

      gplots::heatmap.2(matrix,
                        dendrogram='none',
                        Colv=FALSE,
                        Rowv=TRUE,
                        trace="none",
                        breaks = seq(-0.5, 18.5, 1),
                        col = col_df$col,
                        key = FALSE,
                        cexRow=0.7,
                        add.expr = list(rect(xleft = identity_ranges$start_nt,
                                             xright = identity_ranges$end_nt,
                                             ybottom = par("usr")[3], ytop = par("usr")[4],
                                             border = NA,
                                             col = adjustcolor("blue", alpha = 0.2)),
                                        rect(xleft = aav2_vr_ranges$start_nt,   #VRs
                                             xright = aav2_vr_ranges$end_nt,
                                             ybottom = par("usr")[3], ytop = par("usr")[4],
                                             border = NA, #"red",
                                             density = 20,
                                             col = adjustcolor("red", alpha = 0.5))))

      title(paste0("Variant description of ", library_name),
            line = -2,
            adj = 0.6)
      legend(x="bottomleft",
             legend=rownames(col_df),
             fill=col_df$col,
             title = "AAV serotypes",
             title.adj = 0.2,
             inset=c(-.07, -.07),
             xpd=TRUE,
             box.lwd = 0,
             cex = 0.7)

      dev.off()
      
      graphics.off()
      pdf(file.path(output.dir, "reports", paste0("enriched1_", en1_name), "variant_description_reduced_conserved.pdf"), width=8, height=5)
      # plot.variant.description.conserved(matrix = matrix_nt_reduced, 
      #                                    identity_ranges = as.data.frame(identity_ranges_reduced), 
      #                                    col_df = col_df,
      #                                    library_name = "reduced representative variants")
      
      matrix = matrix_nt_reduced
      identity_ranges = as.data.frame(identity_ranges_reduced)
      col_df = col_df
      library_name = "reduced representative variants\n"
      
      gplots::heatmap.2(matrix, 
                        dendrogram='none', 
                        Colv=FALSE, 
                        Rowv=TRUE, 
                        trace="none", 
                        breaks = seq(-0.5, 18.5, 1),
                        col = col_df$col,
                        key = FALSE, 
                        cexRow=0.7,
                        add.expr = list(rect(xleft = identity_ranges$start_nt, 
                                             xright = identity_ranges$end_nt, 
                                             ybottom = par("usr")[3], ytop = par("usr")[4], 
                                             border = NA, 
                                             col = adjustcolor("blue", alpha = 0.2)),
                                        rect(xleft = aav2_vr_ranges$start_nt,   #VRs
                                             xright = aav2_vr_ranges$end_nt, 
                                             ybottom = par("usr")[3], ytop = par("usr")[4], 
                                             border = NA, #"red", 
                                             density = 20, 
                                             col = adjustcolor("red", alpha = 0.5))))
      
      title(paste0("Variant description of ", library_name), 
            line = -2, 
            adj = 0.6)
      legend(x="bottomleft", 
             legend=rownames(col_df), 
             fill=col_df$col,  
             title = "AAV serotypes", 
             title.adj = 0.2, 
             inset=c(-.07, -.07), 
             xpd=TRUE,
             box.lwd = 0, 
             cex = 0.7)
      
      dev.off()
    }
    
    # if multiple input samples 
    if (length(en1.paths) > 1){
      # DESeq2 normalized counts for between samples comparison
      #########################################################
      #########################################################
      cat("\n\nNormalize counts for between sample comparison with median of rations method from \n")
      cat("Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and \n dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014). \n https://doi.org/10.1186/s13059-014-0550-8")
      cat("\n===================================================================\n\n")
      cat(paste0("Output path: ", file.path(output.dir, "files/counts_normalized.csv"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "files/fold_changes.csv"), "\n"))
      
      counts_normalized <- mor_normalization(data_path = file.path(output.dir, "files/counts.csv"), 
                                             output_path = file.path(output.dir, "files/counts_normalized.csv"))
      fold_changes <- get_fold_change(counts_normalized = counts_normalized, 
                                      output_path = file.path(output.dir, "files/fold_changes.csv"))
      counts_fc <- cbind(counts_normalized, fold_changes)
      
      chimeric_orf <- microseq::readFasta(file.path(output.dir,"files/Chimeric_ORF.fasta"))
      chimeric_orf_names <- as.data.frame(matrix(unlist(stringr::str_split(chimeric_orf$Header, pattern = "/")), 
                                    ncol = 5, byrow = TRUE))[[1]]
      chimeric_orf['X'] <- chimeric_orf_names
      
      counts_mor <- dplyr::inner_join(counts_fc %>%
                                        mutate(Representative = rownames(counts_fc)), 
                                      chimeric_orf[c("Sequence", "X")], by = c("Representative" = "X")) 
      write.csv(counts_mor, file.path(output.dir, "files/counts_mor.csv"), row.names = F, quote = F)
      
      cat("\n\nBetween sample comparison plots \n")
      cat("===================================================================\n\n")
      ### Radar plots
      ################################################
      ################################################
      cat(paste0("Output path: ", file.path(output.dir, "reports/radar_fold_changes.pdf"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "reports/radar_log_fold_changes.pdf"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "reports/radar_normalized_counts.pdf"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "reports/radar_log_normalized_counts.pdf"), "\n"))
      
      
      color_vector <- setNames(RColorBrewer::brewer.pal(length(colnames(counts_normalized)), "Set1"), 
                               gsub("\\..*", "",colnames(counts_normalized))) 
      
      pdf(file.path(output.dir, "reports/radar_fold_changes.pdf"), width=10, height=6)
      plot.radar(data = counts_fc, line_width = 1.5, type = "Fold changes", 
                 min_fc = 1, quantile_prob = 0.8, color_vector = color_vector,
                 title = "Fold changes of the enriched variants (FC>1) \nwith normalized counts > 0.8 quantile")
      dev.off()
      
      pdf(file.path(output.dir, "reports/radar_log_fold_changes.pdf"), width=10, height=6)
      plot.radar(data = counts_fc, line_width = 1.5, type = "Fold changes", log_2 = T,
                 min_fc = 1, quantile_prob = 0.8, color_vector = color_vector,
                 title = "log2 fold changes of the enriched variants (FC>1) \nwith normalized counts > 0.8 quantile")
      dev.off()
      
      pdf(file.path(output.dir, "reports/radar_normalized_counts.pdf"), width=11, height=6)
      plot.radar(data = counts_fc, line_width = 1.5, type = "Normalized counts",
                 min_fc = 1, quantile_prob = 0.8, color_vector = color_vector,
                 title = "Normalized counts of the enriched variants (FC>1) \nwith normalized counts > 0.8 quantile")
      dev.off()

      pdf(file.path(output.dir, "reports/radar_log_normalized_counts.pdf"), width=11, height=6)
      plot.radar(data = counts_fc, line_width = 1.5, type = "Normalized counts", log_2 = T,
                 min_fc = 1, quantile_prob = 0.8, color_vector = color_vector,
                 title = "log2 normalized counts of the enriched variants (FC>1) \nwith normalized counts > 0.8 quantile")
      dev.off()

      ### MA plots
      ################################################
      ################################################
      cat(paste0("Output path: ", file.path(output.dir, "reports/ma_chimeric_normalized_counts.pdf"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "reports/ma_mean_normalized_counts.pdf"), "\n"))
      cat(paste0("Output path: ", file.path(output.dir, "reports/ma_sample_normalized_counts.pdf"), "\n"))
      
      pdf(file.path(output.dir, "reports/ma_chimeric_normalized_counts.pdf"), width=11, height=6)
      plot.ma(counts_normalized = counts_normalized,
              fold_changes = fold_changes,
              color_vector = color_vector, 
              quantile_prob = 0.8, 
              point_size = 2.3,
              title = "log2FC vs normalized counts in chimeric library",
              xlabel = "Normalized counts in chimeric library")
      dev.off()
      
      pdf(file.path(output.dir, "reports/ma_mean_normalized_counts.pdf"), width=11, height=6)
      plot.ma(counts_normalized = counts_normalized,
              fold_changes = fold_changes,
              color_vector = color_vector, 
              quantile_prob = 0.8, 
              point_size = 2.3,
              title = "log2FC vs mean normalized counts",
              xlabel = "Mean of normalized counts")
      dev.off()
      
      pdf(file.path(output.dir, "reports/ma_sample_normalized_counts.pdf"), width=11, height=6)
      plot.ma(counts_normalized = counts_normalized,
              fold_changes = fold_changes,
              color_vector = color_vector, 
              quantile_prob = 0.8, 
              point_size = 2.3,
              title = "log2FC vs normalized counts in samples",
              xlabel = "Normalized counts in samples")
      dev.off()
    }
  }
} else if (identify) {
  cat("\nUse exploreout folder and enriched libs to find abundant variants\n")
  cat("\nNot implemented yet.")
}


output_paths <- na.omit(output_paths)
write.table(output_paths, file.path(output.dir, "config/output_paths"), 
            sep = " ", quote = F, row.names = F, col.names = F)

save.image(rdata_name)

#load(rdata_name)

sink()

