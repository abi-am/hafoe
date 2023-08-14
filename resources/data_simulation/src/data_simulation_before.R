# First specify the packages of interest
packages = c("stringr", "seqinr", "microseq")

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
      install.packages(x, dependencies = TRUE)
      suppressPackageStartupMessages(library(x, character.only = TRUE, quietly = T))
    }
  }
)

# Run this first (before running hafoe) to generate the input files

library(seqinr)
library(stringr)

dir.create("plots", showWarnings = F)
dir.create("files", showWarnings = F)
fileName <- "input_files/AAV_all16_new.clustal_num"
format <- "clustal"

alignment <- seqinr::read.alignment(fileName, format = format, forceToLower = F)

chimseq_num <- 300
perc_ones_chimeric <- 0.8
aav_num <- alignment$nb
fragment_size_min <- 100 
fragment_size_max <- 700 
min_orf <- 600
aln_len <- nchar(alignment$seq[1])


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

biased_aav_nums <- c(rep(2, 25), rep(14, 20), rep(6, 17), rep(7, 10), rep(9, 7), 
                     rep(8, 6), rep(13, 6), rep(1, 5), rep(5, 2), rep(3, 1), 
                     rep(15, 1))

for (d in 1:5){
  print(paste0("Running on data "), d)
  dir.create(paste0("input_files/data", d), showWarnings = F)
  dir.create(paste0("files/data", d), showWarnings = F)
  
  seed <- d
  
  seqs <- c()
  seq_labels <- c()
  
  i <- 1
  while(i <= chimseq_num){
    #for (i in seq(1, chimseq_num)){
    cursor <- 1
    seq <- ""
    seq_label_list <- c()
    
    while (cursor < aln_len) {
      
      #choose random serotype 
      sero_idx <- sample(biased_aav_nums, 1, replace=T)  #uniform
      #choose random cut position
      cut <- sample((cursor + fragment_size_min) : (cursor + fragment_size_max), 1, replace=F) #cursor:aln_len
      
      #if cut position is near the end, take the fragment until the end
      if (aln_len - cut <= fragment_size_min){
        cut <- aln_len
      }
      
      #append seq, to generate the chimeric sequence
      seq <- paste0(seq, substring(alignment$seq[which(alignment$nam == serotype_name[sero_idx])], cursor, cut))  #sero_idx -> 
      #append seq_label, to keep track of the chimeric sequence composition
      seq_label_list <- c(seq_label_list, rep(sero_idx, (cut - cursor + 1)))
      
      #update cursor
      cursor <- cut
      
    }
    
    #remove gap symbols from sequence and label list
    gap_positions <- unlist(gregexpr('-', seq))
    seq <- gsub("-", "", seq)
    
    seq_label_list <- seq_label_list[- gap_positions]
    
    seq_label <- paste(seq_label_list, collapse = " ")
    
    #check for the ORF length
    pos <- ORFik::findORFs(seq, startCodon = c("ATG"), minimumLength = min_orf)
    neg <- ORFik::findORFs(microseq::reverseComplement(seq), startCodon = c("ATG"), minimumLength = min_orf)
    if (length(pos) > 0 || length(neg) > 0){
      # print(pos)
      # print(neg)
      seqs <- c(seqs, seq)
      seq_labels <- c(seq_labels, seq_label)
      i <- i + 1
    } 
  }  
  
  df <- data.frame(chimeric_seq = seqs,
                   composition = seq_labels)
  
  #CHIMERIC
  
  chimeric_library_initial <- df[c('chimeric_seq', 'composition')]
  index <- c(1:nrow(chimeric_library_initial))
  chimeric_library_initial[, "index"] <- index
  chimeric_library_initial[, "X"] = paste0("AAV.", 100000 + as.numeric(chimeric_library_initial$index))
  chimeric_library_initial_fa <- chimeric_library_initial[c("X", "chimeric_seq")]
  colnames(chimeric_library_initial_fa) <- c("Header", "Sequence")
  
  microseq::writeFasta(chimeric_library_initial_fa, paste0("input_files/data", d, "/Chimeric_lib_simulated_initial.fasta"))
  
  
  chimeric_true_labels <- chimeric_library_initial[c('X', 'composition', 'chimeric_seq')] 
  colnames(chimeric_true_labels) <- c('X', 'Composition', 'Sequence')
  
  write.csv(chimeric_true_labels, paste0("input_files/data", d, "/Chimeric_lib_simulated_labels.csv"), row.names = F)
  
  
  # generate counts for chimeric sequences for PacBio sequencing simulation
  #########################################################################
  
  set.seed(seed)
  chimeric_count0 <- rnorm(100000, mean=1, sd=50)
  chimeric_count0 <- chimeric_count0[chimeric_count0 >= 1.5]
  
  set.seed(seed)
  not_one <- round(sample(chimeric_count0, chimseq_num*(1-perc_ones_chimeric)))
  chimeric_count <- rep(1, chimseq_num)
  set.seed(seed)
  not_one_idx <- sample(1:chimseq_num, chimseq_num*(1-perc_ones_chimeric))
  chimeric_count[not_one_idx] <- not_one
  
  # fig <- plotly::plot_ly(x = counts, type = "histogram", xbins = list(size=1)) 
  # fig
  
  write.table(chimeric_count, 
              paste0("files/data", d, "/chimeric_counts.csv"), 
              sep=",", 
              quote = F, 
              row.names = F, 
              col.names=FALSE)
  
  
  #ENRICHED
  ###
  temp <- readFasta(paste0("input_files/data", d, "/Chimeric_lib_simulated_initial.fasta"))
  temp_counts <- read.csv(paste0("files/data", d, "/chimeric_counts.csv"), header = F)
  chimeric_enriched_counts <- data.frame(name = temp$Header,
                                         chimeric_count = temp_counts$V1)
  
  set.seed(seed)
  fraction_change0 <- rnorm(1000, mean = -1, sd = 0.5)
  fraction_change <- fraction_change0[fraction_change0 >= -1]
  set.seed(seed)
  fraction_change <- sample(fraction_change, chimseq_num)
  fraction_change <- round(fraction_change, 2)
  
  enriched_count <- floor((1+fraction_change) * chimeric_enriched_counts$chimeric_count)
  chimeric_enriched_counts['enriched_count'] <- enriched_count
  
  
  write.table(enriched_count, 
              paste0("files/data", d, "/enriched_counts.csv"), 
              sep=",", 
              quote = F, 
              row.names = F, 
              col.names=FALSE)
  
}
