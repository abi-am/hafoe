
#plot sequence length distribution histogram
plot.histogram <- function(seq_data, library_name = "", variant_orf = "", additional_info = "", bin_size = 20){
  fig <- plotly::plot_ly(x = nchar(seq_data), type = "histogram", xbins = list(size=bin_size)) 
  fig <- fig %>% layout(title= list(text = paste0(library_name, " library ", variant_orf, " length distribution", additional_info))
                        #,xaxis = list(range=c(0,1000))
  )
  fig
}

plot.histogram_ <- function(x, library_name = "", variant_orf = "", additional_info = "", bin_size = 20){
  fig <- plotly::plot_ly(x = x, type = "histogram", xbins = list(size=bin_size)) 
  fig <- fig %>% layout(title= list(text = paste0(library_name, " library ", variant_orf, " length distribution", additional_info))
                        #,xaxis = list(range=c(0,1000))
  )
  fig
}


library(dplyr)
deduplicate <-  function(file, out_path) {
  ###
  # function which 
  # 1. takes as an input the fasta/fastq file (file),
  # 2. gets the counts of sequences that are in the fasta file (duplicates),
  # 3. creates the reverse-complements of each sequence that is in the fasta_file and stores it in a new dataframe (duplicates_rc), 
  # 4. combines original sequences' dataframe and reverse-complemented sequences' dataframe (duplicates_both)
  # 5. duplicates_both has columns Sequence, 
  # Count (how many of that sequence are there in the initial fasta file), 
  # Forward (1 if original sequence from fasta file, 0 if reverse-complement of it), 
  # Name (number uniquely identifying original sequences. Original sequences (in duplicates dataframe) have 
  # the same names with their reverse-complements(in duplicates_rc dataframe)) 
  # 6. group by sequence, sum up Count, and take min Name (duplicates_new1), then group by name and take max Forward (duplicates_new2).  
  # 7. returns the final dataframe containing sequences and counts (duplicates_final)
  ###
  
  #table() function shows how many times each sequence occurs in fasta_file
  #convert to dataframe and add column names for convenience 
  duplicates <- as.data.frame(table(file$Sequence))
  colnames(duplicates) <- c("Sequence", "Count")
  duplicates$Sequence <- as.character(duplicates$Sequence)
  #Order by Count of sequence in decreasing order
  duplicates <- duplicates[order(duplicates$Count, decreasing = T),]
  #rewrite rownames 
  rownames(duplicates) <- 1:nrow(duplicates)
  #create new column called Name and which is basically the index of corresponding row 
  #(integer from 1 to number of rows of duplicates dataframe)
  duplicates['Name'] <- 1:nrow(duplicates) 
  #create new column called Forward having value 1 for each row
  duplicates['Forward'] <- 1
  
  #new dataframe containing reverse-complemented sequences of original sequences' dataframe (duplicates), 
  #count is the same as for corresponding sequence in original sequences' dataframe, 
  #Forward is 0 for each row, Name is the same as for corresponding sequence in original sequences' dataframe
  duplicates_rc <- data.frame(Sequence = microseq::reverseComplement(duplicates$Sequence), 
                               Count = duplicates$Count, 
                               Forward = 0, 
                               Name = duplicates$Name)
  
  #combine both into a new dataframe 
  duplicates_both <- rbind(duplicates, duplicates_rc)
  
  
  duplicates_new1 <- duplicates_both %>%
    group_by(Sequence) %>%            #group by Sequence, which means the functions after this will work on the groups of Sequences
    mutate(Count = sum(Count)) %>%    #overwrite the column Count to be the sum of counts of rows having the same sequence
    filter(Name == min(Name)) %>%     #keep only rows with minimum name for each unique sequence
    arrange(desc(Count))              #sort descendingly by Count
  
  duplicates_new2 <- duplicates_new1 %>%
    group_by(Name) %>%                #group by Name
    filter(Forward == max(Forward))   #keep only rows having maximum Forward value which is 1, to keep original sequences
  
  #take only columns Sequence and Count
  duplicates_final <- duplicates_new2[c("Sequence", "Count")]
  index <- c(1:nrow(duplicates_final))
  duplicates_final[, "index"] <- index
  duplicates_final[, "X"] = paste0("AAV.", 100000 + as.numeric(duplicates_final$index))
  
  return(duplicates_final)
  #write into a csv file
  # write.csv(duplicates_final,
  #           out_path,
  #           row.names = FALSE)
}


write.fastq <- function(fastq_file, output_path){
  c = 1
  while (c <= nrow(fastq_file)){
    cat(paste0("@", fastq_file$Header[c]),
        fastq_file$Sequence[c],
        "+",
        fastq_file$Quality[c],
        sep = "\n",
        file = output_path,
        append = TRUE)
    c = c + 1
  } 
}

#filter by minimum ORF length and max sequence length (ORF length >= min_orf (600 aa) and sequence length <= max_seq(3000))
filter.by.orf <- function(file, file_path, min_orf = 600, max_seq = 3000, format = "fasta", library_name = "", output_paths){
  cat(paste0("File: ", file_path, "\n"))
  cat(paste0("Minimum ORF length: ", min_orf, "\n"))
  cat(paste0("Maximum sequence length: ", max_seq, "\n"))
  cat(paste0("Number of sequences before filtering: ", nrow(file), "\n"))
  
  dna <- readDNAStringSet(file_path, format = format)
  dna_rc <- Biostrings::reverseComplement(dna)
  
  #Positive strand
  #minimumLength -> number of amino acids not bases
  pos <- findORFs(dna, longestORF = TRUE, startCodon = c("ATG"), minimumLength = min_orf)
  #Get sequences 
  gr_pos <- unlist(pos, use.names = TRUE)
  gr_pos <- GRanges(seqnames = names(dna)[as.integer(names(gr_pos))], ranges(gr_pos), strand = "+")
  names(gr_pos) <- make.names(gr_pos@ranges@NAMES, unique=TRUE)
  gr_pos_df <- as.data.frame(gr_pos)
  #Give proper names:
  #Ex. m54278_211112_013512/4391649/ccs/ORF/1706-3910/2205/+   
  names(gr_pos) <- paste0(gr_pos_df$seqnames, "/", "ORF/", gr_pos_df$start, "-", gr_pos_df$end, "/", gr_pos_df$width, "/+")
  orf_seq_pos <- Biostrings::getSeq(dna, gr_pos)
  
  orf_info_df_pos <- as.data.frame(gr_pos)
  orf_seq_df_pos <- as.data.frame(orf_seq_pos)
  
  #Negative strand
  neg <- findORFs(dna_rc, longestORF = TRUE, startCodon = c("ATG"), minimumLength = min_orf)
  if (length(neg) > 0){
  #Get sequences
    gr_neg <- unlist(neg, use.names = TRUE)
    gr_neg <- GRanges(seqnames = names(dna_rc)[as.integer(names(gr_neg))], ranges(gr_neg), strand = "+") #strand = "+" to get ORFs instead of ORF_rc
    names(gr_neg) <- make.names(gr_neg@ranges@NAMES, unique=TRUE)
    gr_neg_df <- as.data.frame(gr_neg)
    #Give proper names: 
    names(gr_neg) <- paste0(gr_neg_df$seqnames, "/", "ORF/", gr_neg_df$start, "-", gr_neg_df$end, "/", gr_neg_df$width, "/-")
    orf_seq_neg <- Biostrings::getSeq(dna_rc, gr_neg)
    
    orf_info_df_neg <- as.data.frame(gr_neg)
    orf_seq_df_neg <- as.data.frame(orf_seq_neg)
  } else {
    orf_info_df_neg <- c()
    orf_seq_df_neg <- c()
  }
  
  #Combine pos & neg
  orf_info_df_all <- rbind(orf_info_df_pos, orf_info_df_neg)
  orf_seq_df_all <-  rbind(orf_seq_df_pos, orf_seq_df_neg)
  
  #Take max ORF (group by seqname, take the row with max orf length) (some variants have multiple orfs > 1.8k)
  group <- as.data.table(orf_info_df_all, keep.rownames = TRUE)
  orf_info_df_all_max <- group[group[, .I[width == max(width)], by=seqnames]$V1]
  #Check for duplicates
  #which(duplicated(orf_info_df_all$seqnames))
  
  #Remove duplicate rows by seqname leaving only first: leave only one ORF of same variant having same length
  orf_info_df_all_max <- orf_info_df_all_max[!duplicated(orf_info_df_all_max$seqnames),]
  
  orf_seq_df_all['name'] <- rownames(orf_seq_df_all)
  orf_seq_df_all_max <- orf_seq_df_all[rownames(orf_seq_df_all) %in% orf_info_df_all_max$rn,]
  
  #Contains orfs
  fasta_orf_seq <- orf_seq_df_all_max[c('name', 'x')]
  colnames(fasta_orf_seq) <- c("Header", "Sequence")
  
  #Contains original sequences
  fasta_orf <- dplyr::inner_join(file, orf_info_df_all_max, by = c("Header" = "seqnames")) 
  fasta_orf <- fasta_orf[c('rn', 'Sequence')]  
  colnames(fasta_orf) <- c('Header', 'Sequence')
  
  #Filter by sequence length (max_seq)
  fasta_orf_filtered <- fasta_orf[nchar(fasta_orf$Sequence) <= max_seq, ] 
  fasta_orf_seq_filtered <- fasta_orf_seq[fasta_orf_seq$Header %in% fasta_orf_filtered$Header,]
  
  #Store filtered orf sequences into fasta file
  writeFasta(fasta_orf_seq_filtered, file.path(output.dir, "files", paste0(library_name, "_ORF.fasta")))

  if (format == "fastq"){
    #Add sequence original name, orf start, end positions as columns
    orf_seq_df_all_max$name_original <- sub("/ORF.*", "", orf_seq_df_all_max$name)
    name_splitted <- stringr::str_split(orf_seq_df_all_max$name, pattern = "/", simplify = TRUE) 
    ranges <- apply(name_splitted, 1, FUN=function(x) x[x != ""][length(x[x != ""])-2] )  #[,5]
    start <- stringr::str_split(ranges, pattern = "-", simplify = TRUE)[,1]
    end <- stringr::str_split(ranges, pattern = "-", simplify = TRUE)[,2]
    orf_seq_df_all_max['start'] <- start
    orf_seq_df_all_max['end'] <- end
    
    #Contains orfs
    fastq_orf_seq <- inner_join(file, orf_seq_df_all_max, by = c("Header" = "name_original")) 
    fastq_orf_seq['quality_orf'] <- substring(fastq_orf_seq$Quality, fastq_orf_seq$start, fastq_orf_seq$end)
    fastq_orf_seq <- fastq_orf_seq[c('name', 'x', 'quality_orf')]  
    colnames(fastq_orf_seq) <- c('Header', 'Sequence', 'Quality')
    
    fastq_orf_seq_filtered <- fastq_orf_seq[fastq_orf_seq$Header %in% fasta_orf_filtered$Header,]
    
    # write.fastq(fastq_file = fastq_orf_seq_filtered, 
    #             output_path = file.path(output.dir, "files", paste0(library_name, "_ORF.fastq")))
    
    #quote error
    #writeFastq(fastq_orf_seq_filtered, file.path(output.dir, "files", paste0(library_name, "_ORF.fastq")))
  }
  
  
  cat(paste0("Number of sequences after filtering: ", nrow(fasta_orf_seq_filtered), "\n"))
  cat("ORF length summary statistics\n")
  print(summary(nchar(fasta_orf_seq_filtered$Sequence)))
  cat("\n")
  
  orf_summary <- data.frame(info = c("Minimum ORF length", "Maximum sequence length", "Number of sequences before filtering", "Number of sequences after filtering"),
                                     value = c(min_orf, max_seq, nrow(file), nrow(fasta_orf_seq_filtered)))
  write.csv(orf_summary, file.path(output.dir, "files", paste0(library_name, "_orf_summary.csv")), row.names = F, quote = F)
  
  # store output paths 
  output_paths <- rbind(output_paths, 
                        c(paste0('final_', library_name, '_orf_fasta'), file.path(output.dir, "files", paste0(library_name, "_ORF.fasta"))))
  output_paths <- rbind(output_paths,
                        c(paste0(library_name, "_orf_summary"), file.path(output.dir, "files", paste0(library_name, "_orf_summary.csv"))))
  
  #Store dataframe outputs in list (for plotting histogram)
  out <- list(fasta_orf_filtered, fasta_orf_seq_filtered, output_paths)  
  
  return(out)
}

choose.new.representatives <- function(sizes_path, members_path, chimeric_library, fasta_orf_seq_filtered_chim, fasta_path, csv_path){
  cat("Choosing a new representative ORF sequence for each cluster based on: \n")
  cat("\t1. greatest sequence abundance in the input chimeric library\n")
  cat("\t2. maximum ORF length\n")
  cat(paste0("Output path: ", fasta_path, "\n"))

  #use the generated members_ordered.csv and sizes.csv to separate variants by clusters
  sizes <- read.csv(sizes_path, header = FALSE)
  sizes['sum'] <- cumsum(sizes$V1)
  members <- read.csv(members_path, header = FALSE)

  #for each cluster choose a representative sequence (1.max count, 2.max ORF length)
  #for loop which works for each cluster
  for (i in seq_len(nrow(sizes))){
    #only for 1st cluster
    if (i == 1){
      #make dataframe containing only the variants which are in the 1st cluster (cluster 0)
      #example of a variant name: AAV.211507/ORF/754-2958/2205/-      I split by / to get   AAV.211507  ORF  754-2958  2205  -     this info in separate columns, to be able to use width/length of orfs
      temp <- as.data.frame(stringr::str_split(members[1:sizes$sum[1],], pattern = "/", simplify = TRUE))
      colnames(temp) <- c("name", "ORF", "pos", "width", "strand")


      #join to initial data by variant name to also take the count column
      temp <- dplyr::inner_join(temp, chimeric_library[c('X', 'Count')],  by = c("name" = "X"))
      #take only the row(s)/variant(s) having max count
      max_count <- temp[temp$Count == max(temp$Count),]
      #from them take only the row(s)/variant(s) having max orf length
      max_width <- max_count[max_count$width == max(max_count$width),]
      #make a dataframe rep containing cluster number, representative's name, its orf length and abundance count
      rep <- data.frame(Cluster = i-1,
                        Representative = max_width[1, 'name'], #[1,'name'] 1 means I take only 1st row just in case multiple having max count
                        Width = max_width[1,'width'],
                        Count = max_width[1, 'Count'])
      #make df of cluster members
      cluster_members <- data.frame(Cluster = i-1,
                                    Representative = max_width[1, 'name'])
      cluster_members$Members <- list(temp$name)
    }
    #for the rest of clusters
    else{
      #same thing for the rest of clusters (cluster 1, then cluster 2, etc.)
      temp <- as.data.frame(stringr::str_split(members[(sizes$sum[i-1]+1) : sizes$sum[i], ],
                                            pattern = "/", simplify = TRUE))
      colnames(temp) <- c("name", "ORF", "pos", "width", "strand")

      temp <- dplyr::inner_join(temp, chimeric_library[c('X', 'Count')],  by = c("name" = "X"))
      max_count <- temp[temp$Count == max(temp$Count),]
      max_width <- max_count[max_count$width == max(max_count$width),]
      #populates the dataframe rep created above for the 1st cluster
      rep <- rbind(rep,
                   c(i-1,
                     max_width[1, 'name'],
                     max_width[1, 'width'],
                     max_width[1, 'Count']))
      
      cluster_members <- rbind(cluster_members,
                               c(i-1,
                                 max_width[1, 'name'],
                                 NA))
      cluster_members$Members[i] <- list(temp$name)
      
    }
  }
  
  
  #from initial fasta file take only the representative variants
  fasta_orf_seq_filtered_names <- stringr::str_split(fasta_orf_seq_filtered_chim$Header, pattern = "/", simplify = TRUE)[,1]
  
  rep_fasta <- fasta_orf_seq_filtered_chim[fasta_orf_seq_filtered_names %in% rep$Representative,]
  
  #create a new fasta file containing representatives
  microseq::writeFasta(rep_fasta,
                       fasta_path)
  write.table(rep$Representative,
              csv_path,
              row.names = F,
              col.names = F,
              quote = F,
              sep = ",")
  
  cluster_members$Members_char <- sapply(cluster_members$Members, paste, collapse = " ")
  write.csv(cluster_members[c(-3)], 
            file.path(output.dir, "files/clstr_chimeric_lib/clustering_info.csv"),
            quote = F, 
            row.names = F)
  return(cluster_members)
  
}

#change based on fastq file
make.fastq.files <- function(read_length, fasta.path, step_size, library_name = ""){
  cat("Chopping each representative sequence into overlapping fragments of fixed size and storing into fastq files \n")
  cat(paste0("Fragment length: ", read_length, "\n"))
  cat(paste0("Step size: ", step_size, "\n"))
  cat(paste0("Output path: ", file.path(output.dir, "files/variant_description", library_name, "fastq"), "\n"))
  if (step_size > read_length){
    cat("Warning: step size is greater than fragment length, some regions of the sequence will be skipped!")
  }

  data <- microseq::readFasta(fasta.path)  
  if (str_count(data$Header, '_')[1] == 4){
    data$name <- matrix(unlist(stringr::str_split(data$Header, "/")), ncol = 5, byrow = TRUE)[,1]
  } else {
    data$name <- data$Header
  }  
    
  dir.create(file.path(output.dir, "files/variant_description"), showWarnings = F)
  dir.create(file.path(output.dir, "files/variant_description", library_name), showWarnings = F)
  dir.create(file.path(output.dir, "files/variant_description", library_name, "fastq"), showWarnings = F)
  splitter <- function(rl, ss, c){
    pos = 0     
    read_counter = 10001 #10^nchar(as.character(ceiling(max(nchar(data$Sequence))/ss)))+1
    
    # GENERAL CASE
    while(pos < nchar(data[c, "Sequence"])){  #fragments of length less than rl may be generated at the end of the sequence
      cat(paste0("@", data[c, "Header"], "_", read_counter),
          substring(data[c, "Sequence"], pos + 1, pos + as.numeric(rl)),
          "+",
          strrep("h", nchar(substring(data[c, "Sequence"], pos + 1, pos + as.numeric(rl)))),
          sep = "\n",
          file = file.path(output.dir, "files/variant_description", library_name, "fastq",
                           paste0(data[c,"name"], ".fastq")),
          append = TRUE)
      pos = pos + as.numeric(ss)
      read_counter = read_counter + 1
    }
  }
  
  c = 1
  while (c <= nrow(data)){
    splitter(read_length, step_size, c)
    c = c + 1
  }
}



#parents_df change general case from parents fasta or separate input paramenter from .sh
neighbor.joining <- function(parents_df, path_fasta, path_mapped, 
                             path_unmapped, read_length, step_size,
                             vd_criterion, library_name = ""){
  start_time <- Sys.time()
  
  
  cat(paste0("Output path: ", file.path(output.dir, 
                                        "files/variant_description", 
                                        library_name, 
                                        paste0(library_name, "_variant_description.csv")), "\n"))
  
  #make directory for output file if not created yet
  dir.create(file.path(output.dir, "files/variant_description"), showWarnings = F)
  dir.create(file.path(output.dir, "files/variant_description", library_name), showWarnings = F)
  
  data <- microseq::readFasta(path_fasta)  
  data <- data[order(data$Header),] #order data by Header in alphabetical order
  
  #path_mapped - path of directory containing csv files generated from bam files(only mapped reads). 
  #In these csv files there are 2 columns(read name and serotype name)
  
  #list of file names in that directory
  files_mapped <- list.files(path = path_mapped, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  files_mapped <- sort(files_mapped, decreasing = F) #sort names alphabetically
  
  #path_unmapped - path of directory containing csv files generated from fastq files(only unmapped reads). 
  #In these csv files there is 1 column(read name)
  
  #list of file names in that directory
  files_unmapped <- list.files(path = path_unmapped, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  files_unmapped <- sort(files_unmapped, decreasing = F) #sort names alphabetically
  
  files_names <- sub('\\.csv$', '', list.files(path = path_mapped, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE))
  #binds these two lists of files by column (resulting 2 column and each row corresponding to one variant). Needed to connect mapped and unmapped files for each variant
  file_pairs = cbind(files_mapped, files_unmapped, files_names)  
  
  isEmpty <- function(x) {
    return(length(x)==0)
  }
  
  
  #empty matrix which will be filled in the for loop
  m <- matrix(data=list(),
              nrow=nrow(file_pairs),
              ncol=2)
  read <- c()
  options <- list()
  
  #for loop iterating over all query variants' mapped and unmapped files simultaneously 
  for (j in seq_len(nrow(file_pairs))) {
    #there where some empty files in mapped and unmapped files, so I used if/else statements to ignore them
    #MAPPED
    if (file.size(file.path(path_mapped, file_pairs[[j, 1]])) > 0){
      
      #open csv file, do not use first row as header, use tab to separate columns. read as matrix
      data_mapped <- read.csv(file.path(path_mapped, 
                                        file_pairs[[j, 1]]), 
                                        header = F, 
                                        sep = "\t")
      #specify column names of matrix
      colnames(data_mapped) <- c("variant_read", "serotype", "aln_score") #"aln_score"
      data_mapped[,"aln_score"] <- gsub(".*:", "", data_mapped[,"aln_score"])
      
      #extract read numbers from the first column of data_mapped matrix (ex: AAV.100001.10003 --> 3)
      read_num <- as.numeric(substring(data_mapped[,"variant_read"], 
                                       nchar(data_mapped[,"variant_read"]) - 5 + 1, 
                                       nchar(data_mapped[,"variant_read"]))) - 10000
      
      #add this list of read numbers to data_mapped matrix as a column
      data_mapped <- cbind(data_mapped, read_num)
    } else {
      read_num <- c(-1)
    }
    
    #UNMAPPED
    #if file not empty
    if (file.size(file.path(path_unmapped, file_pairs[[j, 2]])) > 0){ 
      #same for unmapped files
      data_unmapped <- read.csv(file.path(path_unmapped, file_pairs[[j, 2]]), 
                                          header = F, 
                                          sep = "\t")
      
      colnames(data_unmapped) <- c("variant_read")
      
      read_num_un <- as.numeric(substring(data_unmapped[,"variant_read"], 
                                          nchar(data_unmapped[,"variant_read"])- 5 + 1, 
                                          nchar(data_unmapped[,"variant_read"]))) - 10000
    } else {
      read_num_un <- c(-1)
    }
    

    variant_length <- nchar(data[j,]$Sequence) #order for j should be the same in data and file_pairs
    v <- list()
    rn <- list()
    for (p in seq_len(variant_length)){
      rn[p] <- c()
      serotypes_all <- c()
      for (z in 1: max(max(read_num), max(read_num_un))){
        if (p %in% seq(((z-1)*step_size)+1, ((z-1)*step_size)+read_length)){ #read start end positions in variant sequence
          
          #if the number of reads is present in the read_num list (list for mapped reads)
          if (z %in% read_num) {
            #then the vector v is appended by the serotype number in corresponding position 
            s <- data_mapped[data_mapped[,"read_num"] == z, "serotype"]
            s_num <- parents_df[serotype_name %in% s, "serotype_num"]
            
          }
          #if the number of reads is present in the read_num_un list(list for unmapped reads)
          if (z %in% read_num_un) {
            #then  the vector v is appended by 0 in corresponding position 
            s_num = 0 
          } 
          
          serotypes_all <- c(serotypes_all, s_num)
          serotypes_all <- unique(serotypes_all)
          
          #remove 0 if there are serotypes identified for that position
          if (length(serotypes_all) > 1){
            serotypes_all <- serotypes_all[serotypes_all != 0]
          }
          
          v[[p]] <- serotypes_all
          rn[[p]] <- c(unlist(rn[p]), z)
        }
      }
    }
    
    
    #leave only seros with the best aln score for that read
    v_max_aln_score <- list()

    #default min score thresh in bowtie2: -0.6 + -0.6 * L , where L is the read length
    min_score <- round(-0.6 + -0.6 * read_length)
    data_mapped$read_num <- as.numeric(data_mapped$read_num)
    data_mapped$aln_score <- as.numeric(data_mapped$aln_score)
    data_mapped["aln_score_shifted"] <- data_mapped$aln_score + abs(min_score)
    data_mapped <- merge(data_mapped, parents_df, by.x = 'serotype', by.y ='serotype_name')
    
    for (i in seq_len(length(v))){
      
      #if there is at least one mapped read for the i-th position
      if (any(rn[[i]] %in% read_num)) {
        all_reads <- data_mapped[data_mapped$read_num %in% rn[[i]],]
        #all_reads_max <- all_reads[all_reads$aln_score == max(all_reads$aln_score),]
        
        #!!!dplyr package problem, doesn't recognize %>% 
        try(detach("package:plyr", unload = TRUE), silent = T)
        library(dplyr)
        
        score_sum_per_sero <- all_reads %>% 
          group_by(serotype_num) %>%
          summarise(
            score_sum = sum(aln_score_shifted),
            score_avg = mean(aln_score_shifted)
          )
  
        if (vd_criterion == "sum"){
          sero_w_max_score <- score_sum_per_sero[score_sum_per_sero$score_sum == max(score_sum_per_sero$score_sum), ]$serotype_num
        } else if (vd_criterion == "avg"){
          sero_w_max_score <- score_sum_per_sero[score_sum_per_sero$score_avg == max(score_sum_per_sero$score_avg), ]$serotype_num
        }
      #otherwise, if no mapped read for the i-th position
      } else {
        sero_w_max_score <- v[[i]] #which should be 0
      }
      v_max_aln_score[[i]] <- sero_w_max_score
    }
    
    v <- v_max_aln_score
    
    # get unique consecutive values from v and lengths of them in a separate vector
    seen <- v[[1]]
    v_consecutive_unique <- v[1]
    count <- 1 
    consecutive_counts <- c()
    for (i in v[2:length(v)]){
      if (!identical(i, seen)){
        v_consecutive_unique <- append(v_consecutive_unique, list(i))
        consecutive_counts <- c(consecutive_counts, count)
        seen <- i
        count <- 1
      } else {
        count <- count + 1 
      }
    }
    consecutive_counts <- c(consecutive_counts, count)
    

    ###neighbor-aware serotype identification
    v <- v_consecutive_unique
    
    # check only the positions having multiple serotypes that need to be resolved
    multiple_idx <- which(lapply(v, length) > 1)

    for (i in multiple_idx){
      #if there is intersect and they are not the same assign the intersect to both neighbors
      #compare with the left neighbor
      left_intersect <- c()
      right_intersect <- c()
      
      if ((i != 1) & (!setequal(v[i], v[i-1]))){
        if(!isEmpty(intersect(v[[i]], v[[i-1]]))){
          left_intersect <- intersect(v[[i]], v[[i-1]])
        }
      } else if ((i != 1) & (setequal(v[i], v[i-1]))){
        left_intersect <- v[[i-1]]
      } else if (i == 1){
        left_intersect <- c()
      }
      #then compare with the right neighbor
      if ((i != length(v)) & (!setequal(v[i], v[i+1]))){
        if (!isEmpty(intersect(v[[i]], v[[i+1]]))){
          right_intersect <- intersect(v[[i]], v[[i+1]])
        } 
      } else if ((i != 1) & (setequal(v[i], v[i+1]))){
        right_intersect <- v[[i+1]]
      } else if (i == length(v)){
        right_intersect <- c()
      }
      
      # change only the i-th position, not the neighbors
      if (!is.null(union(left_intersect, right_intersect))){ #when both left and right intersects are null, v[[i]] shouldn't change
        v[[i]] <- union(left_intersect, right_intersect)
      }
      
    }

    # make v of sequence length
    v <- rep(v, consecutive_counts)
    
    #assign 17 to positions remaining with multiple serotypes
    for (i in 1:length(v)){
      if (length(v[[i]]) > 1){
        v[[i]] <- 17
      }
    }

    #after the vector v is made, the empty matrix mat defined above is appended. 1st column is the name of variant,
    #2nd column is a string of serotype numbers separated by spaces
    m[[j, 1]] <- file_pairs[[j, 3]]
    m[[j, 2]] <- paste(v, collapse = " ")
    
  }
  
  
  end_time <- Sys.time()
  #end_time - start_time
  
  
  s <- stringr::str_split(unlist(m[,2], 1), " ")
  col_num <- max(nchar(data$Sequence))
  for (i in seq_len(length(s))){
    if (length(s[[i]]) < col_num){ 
      s[[i]] <- c(s[[i]], rep("18", col_num - length(s[[i]]))) #gap
    }
  }
  
  m_new <- matrix(as.numeric(unlist(s)), ncol = col_num, byrow = TRUE)
  rownames(m_new) <- m[,1]
  
  write.table(m_new, 
              file.path(output.dir, 
                        "files/variant_description", 
                        library_name, 
                        paste0(library_name, "_variant_description.csv")), 
              #sep = ",",
              quote = F, 
              col.names = F)
}



plot.cluster.size.distribution <- function(cluster_members, size_thresh, similarity_thresh, library_name = ""){
  cluster_members["Size"] <- unlist(lapply(cluster_members$Members, length))
  data <- cluster_members[-3]
  data <- data[data$Size > size_thresh,]
  
  p <- ggplot2::ggplot(data=data, 
                      aes(x=factor(Representative, data[order(Size, decreasing = T), "Representative"]), 
                          y=Size, 
                          order=Size)) +
        geom_bar(stat="identity")+
        geom_text(aes(label=Size), vjust=-0.3, size=1.5) +
        labs(title=paste0("Cluster size distribution in ", library_name, "\n(Cluster size > ", size_thresh, ", Similarity threshold = ", similarity_thresh,"%)"),
             x ="", y = "Cluster size") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) 
  print(p)
}



plot.variant.description <- function(matrix, col_df, library_name = ""){
  gplots::heatmap.2(matrix, 
                    dendrogram='none', 
                    Colv=FALSE, 
                    Rowv=FALSE, #TRUE 
                    trace="none", 
                    breaks = seq(-0.5, 18.5, 1), 
                    col = col_df$col, 
                    key = FALSE, 
                    cexRow=0.7,
                    margins = c(5,7))
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

}


plot.variant.description.conserved <- function(matrix, col_df, identity_ranges, library_name = ""){
  gplots::heatmap.2(matrix, 
                    dendrogram='none', 
                    Colv=FALSE, 
                    Rowv=TRUE, 
                    trace="none", 
                    breaks = seq(-0.5, 18.5, 1),
                    col = col_df$col,
                    key = FALSE, 
                    cexRow=0.7,
                    margins = c(5,7),
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
}



get.frequency.table <- function(chimeric_library, nj_matrix, cluster_members, output_path){
  #combined counts of chimeric lib representatives
  counts_chim <- data.frame(Representative = NA, 
                            Chimeric.Count = NA)
  
  for (rep in cluster_members$Representative){
    s <- sum(chimeric_library[chimeric_library$X %in% unlist(cluster_members$Members[cluster_members$Representative == rep]), c('Count')])
    counts_chim <- rbind(counts_chim, c(rep, s))
  }
  
  counts_chim <- na.omit(counts_chim)
  counts_chim$Chimeric.Count <- as.numeric(counts_chim$Chimeric.Count)
  
  counts_chim['Chimeric.Percentage'] <- round((100 * (counts_chim$Chimeric.Count / sum(counts_chim$Chimeric.Count))), 2)
  counts_chim <- counts_chim %>% arrange(desc(Chimeric.Count))
  
  write.csv(counts_chim, file.path(output_path, "chimeric_lib_rep_counts.csv"), row.names = FALSE, quote = F)
  
  #serotype frequencies per representative variant
  frequencies = matrix(nrow = nrow(nj_matrix), ncol = 0)
  for (i in (0:18)){
    frequencies = cbind(frequencies, rowSums(nj_matrix == i))
  }
  
  #serotype freq * representative combined count
  frequencies_with_count = matrix(nrow = 0, ncol = ncol(frequencies))
  for (i in rownames(frequencies)){
    frequencies_with_count = rbind(frequencies_with_count, 
                                   counts_chim[counts_chim$Representative == i, "Chimeric.Count"] * frequencies[i,])
  }
  
  frequencies_with_count_df = as.data.frame(frequencies_with_count)
  colnames(frequencies_with_count_df) <- c("no alignment","AAV1","AAV2","AAV3","AAV4","AAV5","AAV6","AAV7","AAV8","AAV9",
                                          "AAV10","AAV11","AAV12","AAV13","AAVrh8","AAVrh10","AAVrh32", "multiple alignment", "gap")
  rownames(frequencies_with_count_df) <- rownames(frequencies)
  
  frequencies_with_count_df <- subset(frequencies_with_count_df, select = -c(gap))
  
  write.csv(frequencies_with_count_df, file.path(output_path, "chimeric_lib_serotype_counts.csv"), quote = F)
  
  #for all representatives
  frequencies_final = colSums(frequencies_with_count_df)
  
  #summary frequency table
  cat("\nAbundance of AAV serotypes in the chimeric library\n")
  serotypes_freq <- as.data.frame(frequencies_final)
  colnames(serotypes_freq) <- c("Freq")
  serotypes_freq['Freq(%)'] <- round(serotypes_freq$Freq*100/sum(serotypes_freq$Freq), 2)
  serotypes_freq <- serotypes_freq[order(serotypes_freq$`Freq(%)`, decreasing = T), ]
  print(serotypes_freq[-1])
  
  write.csv(serotypes_freq, file.path(output_path, "chimeric_lib_serotype_frequency.csv"), quote = F)
  return(serotypes_freq)
}


plot.serotype.frequency <- function(serotypes_freq, col_df, library_name = ""){
  serotypes_freq['Name'] <- rownames(serotypes_freq)
  col_ordered <- col_df[serotypes_freq[order(serotypes_freq$`Freq(%)`, decreasing = T), "Name"],]
  
  p <- ggplot2::ggplot(data=serotypes_freq, 
                        aes(x=factor(Name, serotypes_freq[order(`Freq(%)`, decreasing = T), "Name"]), 
                            y=`Freq(%)`, 
                            order=`Freq(%)`)) +
          geom_bar(stat="identity", fill=col_ordered)+
          geom_text(aes(label=`Freq(%)`), vjust=-0.3, size=3.5) + 
          labs(title=paste0("Distribution of AAV serotypes in ", library_name),
               x ="", y = "Frequency (%)") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  print(p)
}
  

plot.cluster.abundance <- function(file_path, size_thresh, similarity_thresh, library_name = ""){
  data = read.csv(file_path)
  data <- data[data$Chimeric.Count > size_thresh,]
  
  p <- ggplot2::ggplot(data=data, 
                        aes(x=factor(Representative, data[order(Chimeric.Count, decreasing = T), "Representative"]), 
                            y=Chimeric.Count, 
                            order=Chimeric.Count)) +
          geom_bar(stat="identity")+
          geom_text(aes(label=Chimeric.Count), vjust=-0.3, size=1.5) +
          labs(title=paste0("Cluster abundance in ", library_name, "\n(Combined count > ", size_thresh, ", Similarity threshold = ", similarity_thresh, "%)"),
               x ="", y = "Combined count of variants") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) 
  
  print(p)
}

get.counts.table <- function(chim_file_path, clstr_enriched1_dir, en1_name, output_path){
  if (file.exists(file.path(output_path, "counts.csv"))){
    counts_chim <- read.csv(file.path(output_path, "counts.csv"))
    counts_chim <- subset(counts_chim, select = -c(Sequence))
  } else {
    counts_chim <- read.csv(chim_file_path)
  }
  
  sizes_en1 <- read.csv(file.path(clstr_enriched1_dir, "cluster_sizes.csv"), header = FALSE)
  sizes_en1$V1 <- sizes_en1$V1 - 1
  sizes_en1['sum'] <- cumsum(sizes_en1$V1)
  members <- read.csv(file.path(clstr_enriched1_dir, "members_ordered.csv"), header = FALSE)
  reps <- read.csv(file.path(clstr_enriched1_dir, "representatives.csv"), header = FALSE)
  
  
  rep_size_en1 <- data.frame(Representative = reps$V1, Size = sizes_en1$V1)
  rep_size_en1_ordered <- rep_size_en1[order(rep_size_en1$Size, decreasing = T),]
  row.names(rep_size_en1_ordered) <- NULL
  
  counts_en1 <- data.frame(Representative = stringr::str_split(rep_size_en1_ordered$Representative, pattern = "/", simplify = TRUE)[,1], 
                           Enriched1.Count = rep_size_en1_ordered$Size, 
                           Enriched1.Percentage = round((100 * (rep_size_en1_ordered$Size / sum(rep_size_en1_ordered$Size))), 2))
  
  names(counts_en1)[names(counts_en1) == "Enriched1.Count"] <- paste0("Enriched1_", en1_name, ".Count")
  names(counts_en1)[names(counts_en1) == "Enriched1.Percentage"] <- paste0("Enriched1_", en1_name, ".Percentage")
  
  counts_en1 <- counts_en1[order(counts_en1[, paste0("Enriched1_", en1_name, ".Count")], decreasing = TRUE),]
  
  counts_all <- plyr::join_all(list(counts_chim, counts_en1), by='Representative', type='full')
  counts_all <- counts_all %>% arrange(desc(Chimeric.Count))
  
  #ratio/ enrichment rate
  counts_all['Ratio_En1Chim'] <- round(counts_all[, paste0("Enriched1_", en1_name, ".Percentage")] / counts_all$Chimeric.Percentage, 2)
  counts_all$Ratio_En1Chim[which(!is.finite(counts_all$Ratio_En1Chim))] <- 0
  names(counts_all)[names(counts_all) == "Ratio_En1Chim"] <- paste0("Ratio_En1", en1_name, "_Chim")
  
  
  chimeric_orf <- microseq::readFasta(file.path(output_path,"Chimeric_ORF.fasta"))
  names <- as.data.frame(matrix(unlist(stringr::str_split(chimeric_orf$Header, pattern = "/")), 
                                ncol = 5, byrow = TRUE))[[1]]
  chimeric_orf['X'] <- names
  
  counts_all <- dplyr::left_join(counts_all, chimeric_orf[c("X", "Sequence")], by = c("Representative" = "X"))

  write.csv(counts_all, file.path(output_path, "counts.csv"), quote = F, row.names = F)
}


plot.top.reps.in.enrichedlib.ggplot <- function(counts_file_path, topn_thresh, en1_name, library_name = ""){
  counts <- read.csv(counts_file_path)
  counts <- counts %>% arrange(desc(counts[paste0("Enriched1_", en1_name, ".Count")])) 
  counts <- counts[1:topn_thresh,]
  level_order <- counts[, 'Representative']
  
  p <- ggplot2::ggplot(data=counts, 
                       aes(x=counts[1:topn_thresh, paste0("Enriched1_", en1_name, ".Count")], 
                           y=factor(counts[1:topn_thresh, "Representative"], levels = level_order), 
                           order=counts[,paste0("Enriched1_", en1_name, ".Count")])) +
    geom_bar(stat="identity")+
    geom_text(aes(label=counts[,paste0("Enriched1_", en1_name, ".Count")]), hjust=-0.2, size=2.5) +
    labs(title=paste0("Top ", topn_thresh, " representative sequences after clustering \n", library_name), 
         x = 'Number of reads in the cluster', 
         y = 'Representative sequence') +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) #+
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) 
  print(p)
}



filter.enriched.reduced.reps <- function(counts_file_path, output_path, 
                                         enriched_chim_thresh = 0, enriched_en1_thresh, enriched_ratio_thresh, 
                                         reduced_chim_thresh, reduced_en1_thresh = 0, reduced_ratio_thresh,
                                         en1_name){
  counts <- read.csv(counts_file_path)
  
  enriched <- counts[counts$Chimeric.Percentage >= enriched_chim_thresh & counts[, paste0("Enriched1_", en1_name, ".Percentage")] >= enriched_en1_thresh & counts[, paste0("Ratio_En1", en1_name, "_Chim")] >= enriched_ratio_thresh,]
  enriched <- enriched[order(enriched[, paste0("Ratio_En1", en1_name, "_Chim")], decreasing = T),]
  rownames(enriched) <- NULL
  
  reduced <- counts[counts$Chimeric.Percentage >= reduced_chim_thresh & counts[, paste0("Enriched1_", en1_name, ".Percentage")] >= reduced_en1_thresh & counts[, paste0("Ratio_En1", en1_name, "_Chim")] < reduced_ratio_thresh,]
  #reduced <- na.omit(reduced) 
  reduced <- reduced[order(reduced[, paste0("Ratio_En1", en1_name, "_Chim")], decreasing = T),]
  rownames(reduced) <- NULL
  
  write.csv(enriched, file.path(output_path, "enriched_representatives.csv"), quote = F, row.names = F)
  write.csv(reduced, file.path(output_path, "reduced_representatives.csv"), quote = F, row.names = F)
  
  enriched_fa <- enriched[c("Representative", "Sequence")]
  # limit max 20 sequences
  enriched_fa <- enriched_fa[enriched_fa$Representative %in% enriched[1:min(nrow(enriched), 20),]$Representative,]
  colnames(enriched_fa) <- c("Header", "Sequence")
  microseq::writeFasta(enriched_fa, file.path(output_path, "enriched_representatives.fasta"))
  
  reduced_fa <- reduced[c("Representative", "Sequence")]
  # limit max 20 sequences
  reduced_fa <- reduced_fa[reduced_fa$Representative %in% reduced[1:min(nrow(reduced), 20),]$Representative,]
  colnames(reduced_fa) <- c("Header", "Sequence")
  microseq::writeFasta(reduced_fa, file.path(output_path, "reduced_representatives.fasta"))
  
  out <- list(enriched, reduced)
  return(out)
}

plot.enrichment.tiles <- function(enriched, reduced, type = "Enriched", axis_text_size = 6, en1_name){
  if (type == "Enriched"){
    counts_filtered = enriched
    counts_filtered_other = reduced
  } else if (type == "Reduced"){
    counts_filtered = reduced
    counts_filtered_other = enriched
  } else {
    cat("\nPlease provide valid type: enriched or reduced.")
  }
  
  level_order <- counts_filtered[order(counts_filtered[, paste0("Ratio_En1", en1_name, "_Chim")], decreasing = F), 'Representative']
  mid = (round(min(min(counts_filtered$Chimeric.Percentage), min(counts_filtered_other$Chimeric.Percentage),
                   min(counts_filtered[, paste0("Enriched1_", en1_name, ".Percentage")]), min(counts_filtered_other[, paste0("Enriched1_", en1_name, ".Percentage")])) +
                 max(max(counts_filtered$Chimeric.Percentage), max(counts_filtered_other$Chimeric.Percentage),
                     max(counts_filtered[, paste0("Enriched1_", en1_name, ".Percentage")]), max(counts_filtered_other[, paste0("Enriched1_", en1_name, ".Percentage")])))+1)/2
  
  ggp1 <- ggplot(counts_filtered, aes(factor(Representative, level = level_order), "")) +    
    geom_tile(aes(fill = Chimeric.Percentage, height = 1)) +
    geom_text(aes(label = Chimeric.Percentage)) + 
    scale_fill_gradient2(low = "#1b98e0", mid = "white", high = "#EE4B2B", 
                         midpoint = mid,
                         name = "Chim[%]") +
    ylab("") +
    xlab("") +  
    ggtitle(paste0(type, " representative sequences after clustering \nRatio = Enriched1_", en1_name, "/Chimeric")) + 
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text=element_text(size=axis_text_size), plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5.5))
  
  
  ggp2 <- ggplot(counts_filtered, aes(factor(Representative, level = level_order), "")) +    
    geom_tile(aes(fill = counts_filtered[, paste0("Enriched1_", en1_name, ".Percentage")], height = 1)) +
    geom_text(aes(label = counts_filtered[, paste0("Enriched1_", en1_name, ".Percentage")])) + 
    scale_fill_gradient2(low = "#1b98e0", mid = "white", high = "#EE4B2B", 
                         midpoint = mid,
                         name = "En1[%]") +
    ylab("") +
    xlab("") +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text=element_text(size=axis_text_size), plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5.5))
  
  
  ggp3 <- ggplot(counts_filtered, aes(factor(Representative, level = level_order), "")) +    
    geom_tile(aes(fill = counts_filtered[, paste0("Ratio_En1", en1_name, "_Chim")], height = 1)) +
    geom_text(aes(label = counts_filtered[, paste0("Ratio_En1", en1_name, "_Chim")])) + 
    scale_fill_gradient2(low = "#1b98e0", mid = "white", high = "#EE4B2B", 
                         midpoint = (min(min(counts_filtered[, paste0("Ratio_En1", en1_name, "_Chim")]), min(counts_filtered_other[, paste0("Ratio_En1", en1_name, "_Chim")])) + 
                                       max(max(counts_filtered[, paste0("Ratio_En1", en1_name, "_Chim")]), max(counts_filtered_other[, paste0("Ratio_En1", en1_name, "_Chim")])))/2, 
                         name = "Ratio") +
    ylab("") +
    xlab("") +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text=element_text(size=axis_text_size), plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5.5))
  
  

  print(cowplot::plot_grid(ggp1, ggp2, ggp3, align = "hv", ncol = 1))
  
}

# helper functions called in get.conserved.positions()
# function which finds the ranges given list of numbers
findRanges_helper <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x) == 1) as.character(x) else paste0(x[1], "-", x[length(x)]) #length(x) == 1:2
  }), use.names=FALSE)
}

#duplicate single positions
dup_helper <- function(x) {
  if(length(x) == 1) {
    x = rep(x, 2)
  } else {
    x = x
  }
}

get.conserved.positions <- function(aln_file_path, format = "clustal"){
  alignment <- seqinr::read.alignment(aln_file_path, format = format, forceToLower = F)
  alignment$seq <- gsub("[\r\n\t]", "", alignment$seq)
  
  conn <- file(aln_file_path, open="r")
  linn <- readLines(conn)
  close(conn)
  
  #remove not needed lines
  linn <- linn[which(!(startsWith(linn, "CLUSTAL") | linn == ""))]
  
  #get identity pattern line
  identity_pattern = ""
  name_plus_6 <- max(nchar(alignment$nam)) + 6     #name with max length + 6 spaces
  for (i in seq(alignment$nb + 1, length(linn), alignment$nb + 1)){
    identity_pattern = paste0(identity_pattern, substring(linn[i], name_plus_6 + 1, nchar(linn[i])))
  }
  
  #find positions of * : . these symbols in identity_pattern string
  identity_symbol_idx  <- lapply(strsplit(identity_pattern, ''), 
                                 function(x) which(x == '*' | x == ":" | x == "."))
  identity_symbol_idx <- unlist(identity_symbol_idx)
  
  
  #get the ranges where there are identity symbols
  identity_ranges <- findRanges_helper(identity_symbol_idx)
  
  #remove single positions
  #identity_ranges = identity_ranges[unlist(lapply(identity_ranges, function(x) grepl("-", x)))]
  
  #split by - to separate start and end positions
  identity_ranges_splitted <- str_split(identity_ranges, pattern = "-")
  

  identity_ranges_splitted <- lapply(identity_ranges_splitted, dup_helper)
  
  #make dataframe of identity ranges
  identity_ranges_df = as.data.frame(matrix(unlist(identity_ranges_splitted), 
                                            ncol = 2, byrow = TRUE))
  colnames(identity_ranges_df) = c("start_nt", "end_nt")
  identity_ranges_df$start_nt = as.integer(identity_ranges_df$start_nt)
  identity_ranges_df$end_nt = as.integer(identity_ranges_df$end_nt)
  
  #filter to take only ranges with length >= 10
  #identity_ranges_df = identity_ranges_df[(identity_ranges_df$end_aa - identity_ranges_df$start_aa) >= 10,]
  
  #convert aa positions to nt
  # identity_ranges_df['start_nt'] <- identity_ranges_df$start_aa * 3 - 2
  # identity_ranges_df['end_nt'] <- identity_ranges_df$end_aa * 3 
  
  #for heatmap
  # identity_ranges_df['start_nt_'] <- identity_ranges_df$start_nt / 100 #read_length
  # identity_ranges_df['end_nt_'] <- identity_ranges_df$end_nt /100
  
  out = list(identity_ranges_df, alignment)
  return(out)
}


add.gap.info <- function(alignment, nj_matrix, step_size, nj_by_nt = F){
  aln_length <- nchar(alignment$seq[1])
  #split to characters
  aln_seq <- lapply(alignment$seq , function(y) unlist(strsplit(y, '')))
  aln_seq <- lapply(aln_seq, unlist)
  
  #gap positions in alignment sequences
  aln_gap_positions  <- lapply(aln_seq, 
                               function(x) which(x == '-'))
  
  aln_all_positions <- seq_len(aln_length) 
  
  aln_not_gap_positions <- lapply(aln_gap_positions, 
                                  function(x) setdiff(aln_all_positions, x))
  
  # adds gap positions for one variant
  fun <- function(gap_i, not_gap_i, name_i){
    l <- c()
    l[gap_i] <- 18
    
    if (nj_by_nt){
      for (i in seq_len(length(nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]))) {
        ng_idx <- not_gap_i[i]
        l[ng_idx] <- nj_matrix[name_i, i]
      }
    } else {
      j = 0
      for (i in seq_len(length(nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]))){ 
        ng_idx <- not_gap_i[(j + 1):min(length(not_gap_i), (j+step_size))]   #read_length for no overlap
        l[ng_idx] <- nj_matrix[name_i, i]
        j = j + step_size
      }
    }
    return(l)
  }
  
  list_all_nt <- mapply(fun, gap_i = aln_gap_positions, not_gap_i = aln_not_gap_positions, name_i = alignment$nam, SIMPLIFY = F)
  
  matrix_nt <- matrix(as.numeric(unlist(list_all_nt)), nrow = alignment$nb, byrow = TRUE)
  rownames(matrix_nt) <- alignment$nam
  
  return(matrix_nt)
}

#VP positions from https://www.ncbi.nlm.nih.gov/nuccore/110645916
plot.vp.positions <- function(){
  p <- ggplot() + 
        geom_segment(aes(x=1,xend=2208,y=2.5,yend=2.5), color = "blue", size=1.5) + #	2,203..4,410 Length:	2208 nt
        geom_segment(aes(x=412,xend=2208,y=2,yend=2), color = "green", size=1.5) + #2,614..4,410 Length:	1,797 nt
        geom_segment(aes(x=607,xend=2208,y=1.5,yend=1.5), color = "red", size=1.5) + #2,809..4,410 Length:	1,602 nt
        geom_hline(yintercept = 1, color = "white") + 
        geom_hline(yintercept = 3, color = "white") + 
        annotate("text", x = -61, y = 2.5, label = "VP1") +
        annotate("text", x = 352, y = 2, label = "VP2") +
        annotate("text", x = 547, y = 1.5, label = "VP3") +
        theme_classic() +
        xlab("") + ylab("") + 
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.line.y=element_blank()
        )
  print(p)
}


get.reps.nj.matrix.nt <- function(output_path, path_fasta, nj_matrix, nj_by_nt = F){
  data <- microseq::readFasta(path_fasta)
  ranges <- stringr::str_split(data$Header, pattern = "/", simplify = TRUE)[,3]  
  data['Start_orf'] <- stringr::str_split(ranges, pattern = "-", simplify = TRUE)[,1]
  data['End_orf'] <- stringr::str_split(ranges, pattern = "-", simplify = TRUE)[,2]
  data$Header <- sub("/ORF.*", "", data$Header)
  col_num <- max(nchar(data$Sequence))
  
  if(!nj_by_nt){
    fun <- function(name_i){
      positions <- seq_len(nchar(data[sub("/ORF.*", "", data$Header) == name_i,]$Sequence))
      l <- c()
      j = 0
      for (i in seq_len(length(nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]))){
        idx <- positions[(j + 1):min(length(positions), (j+step_size))]
        l[idx] <- nj_matrix[name_i, i]
        j = j + step_size
      }
      l <- c(l, rep(18, col_num - length(l)))
      return(l)
    }
    
    list_all_nt <- lapply(rownames(nj_matrix), fun)
    reps_nj_mat_nt <- matrix(as.numeric(unlist(list_all_nt)), nrow = nrow(nj_matrix), byrow = TRUE)
    rownames(reps_nj_mat_nt) <- rownames(nj_matrix)
    
    seq_labels <- c()
    for (i in list_all_nt){
      i <- i[i != 18]
      seq_label <- paste(i, collapse = " ")
      seq_labels <- c(seq_labels, seq_label)
    }
  } else {
    seq_labels <- c()
    for (name_i in rownames(nj_matrix)){ 
      temp <- nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]
      seq_label <- paste(temp, collapse = " ")
      seq_labels <- c(seq_labels, seq_label)
    }
  }
  
  chimeric_rep_predicted_labels <- data.frame(X = rownames(nj_matrix), 
                                              Composition = seq_labels)
  chimeric_rep_predicted_labels <- dplyr::inner_join(chimeric_rep_predicted_labels, 
                                                     data[c("Header", "Sequence", "Start_orf", "End_orf")], 
                                                     by = c("X" = "Header")) 
  
  write.csv(chimeric_rep_predicted_labels, file.path(output_path, "Chimeric_rep_predicted_labels.csv"), 
            row.names = F, quote = F)
  
}

mor_normalization <- function(data_path, output_path){
  # function from https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html
  library(dplyr)
  library(tibble)
  
  # prepare input data
  data <- read.csv(data_path)
  rownames(data) <- data$Representative
  data <- data[grepl(".Count", colnames(data))]
  
  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>% 
    rownames_to_column('gene') %>% 
    mutate (gene_averages = rowMeans(log_data)) %>% 
    filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  
  write.csv(manually_normalized, output_path, row.names = T, quote = F)
  return(manually_normalized)
}

get_fold_change <- function(counts_normalized, output_path){
  ref <- colnames(counts_normalized)[grepl("Chimeric", colnames(counts_normalized))]
  not_ref <- colnames(counts_normalized)[colnames(counts_normalized) != ref]

  fc <- data.frame()
  for (i in 1:length(not_ref)){
    fold_change <- round(counts_normalized[not_ref[i]] / counts_normalized[ref], 2)
    colnames(fold_change) <- paste0(gsub("\\..*","",not_ref[i]), ".FC")
    if(ncol(fc) == 0){
      fc <- fold_change
    } else {
      fc <- cbind(fc, fold_change)
    }
  }
  
  write.csv(fc, output_path, quote = F, row.names = T)
  return(fc)
}


plot.radar <- function(data, 
                       line_width = 1.2, 
                       type = "Fold changes",
                       log_2 = FALSE, 
                       min_fc = 1,
                       #min_expression = 10,
                       quantile_prob = 0.8,
                       color_vector,
                       title = ""){

  filtered_variants <- c()
  for (i in gsub("\\..*","",colnames(data))){
    if (!grepl("Chimeric", i)){
      min_expression <- quantile(data[,paste0(i, ".Count")], quantile_prob)[[1]] 
      filtered_variants <- c(filtered_variants,
                             rownames(data[(data[paste0(i, ".Count")] >= min_expression) & (data[paste0(i, ".FC")] >= min_fc),]))
    }
  }
  data <- data[unique(filtered_variants),]
  
  if (type == "Normalized counts"){
    count_cols <- colnames(data)[grepl("Count", colnames(data))] 
    data <- data[count_cols]
    if (log_2){
      data <- log2(data + min(data[data!=0])/2) 
    }
  } else if (type == "Fold changes"){
    fc_cols <- colnames(data)[grepl("FC", colnames(data))]
    data <- data[fc_cols]
    if (log_2){
      data <- log2(data + min(data[data!=0])/2) 
    }
  }
  data <- t(data)
  rownames(data) <-  gsub("\\..*","",rownames(data))
  
  # To use the fmsb package, add 2 lines to the dataframe: the max and min of each variable to show on the plot!
  data <- rbind(rep(round(max(data)), ncol(data)), 
                rep(floor(min(data)),ncol(data)) , data)
  data <- as.data.frame(data)

  color_vector <- color_vector[rownames(data)[!grepl("X",rownames(data))]]
  # plot with default options:
  fmsb::radarchart(data , axistype=1 , 
                   #custom polygon
                   pcol=color_vector, 
                   plwd=line_width, plty=5, #2,5
                   pty=32,
                   #custom the grid
                   cglcol="grey", cglty=3, 
                   axislabcol="grey", cglwd=0.8, 
                   caxislabels=round(seq(floor(min(data)),max(data),max(data)/4)), 
                   seg = length(round(seq(floor(min(data)),max(data),max(data)/4)))-1,
                   calcex=0.9,
                   #custom labels
                   vlcex=0.6,
                   title = title)
  
  # Add a legend
  legend(x = 1.5, y=0.5,  
         legend = gsub("\\..*","",rownames(data[-c(1,2),])), 
         bty = "n", pch=20 , col=color_vector, 
         cex=0.8, pt.cex=3)
  
}


plot.ma <- function(counts_normalized,
                    fold_changes,
                    color_vector, 
                    #thresh = 20,
                    quantile_prob = 0.8,
                    hjust_ = -0.3, 
                    xlabel = "Normalized counts in chimeric library", 
                    title = "", 
                    point_size=2){
  
  counts_normalized['Representative'] <- rownames(counts_normalized)
  counts_normalized_melted <- reshape2::melt(counts_normalized, id=c("Representative"))
  counts_normalized_melted$variable <- sub("\\..*", "",counts_normalized_melted$variable)
  
  fold_changes['Representative'] <- rownames(fold_changes)
  fold_changes_melted <- reshape2::melt(fold_changes, id=c("Representative"))
  fold_changes_melted$variable <- gsub("\\..*", "",fold_changes_melted$variable)
  
  #chimeric
  counts_normalized_chimeric <- counts_normalized[c("Chimeric.Count", "Representative")] 
  colnames(counts_normalized_chimeric) <- c("value", "Representative")
  
  #means of all samples
  counts_normalized_means <- as.data.frame(rowMeans(subset(counts_normalized, select = -c(Representative))))
  colnames(counts_normalized_means) <- c("value")
  counts_normalized_means['Representative'] <- rownames(counts_normalized_means)
  
  if (xlabel == "Normalized counts in chimeric library"){
    df <- merge(counts_normalized_chimeric, fold_changes_melted, by=c('Representative'))
  } else if (xlabel == "Normalized counts in samples"){
    df <- merge(counts_normalized_melted, fold_changes_melted, by=c('Representative', 'variable'))
  } else if (xlabel == "Mean of normalized counts"){
    df <- merge(counts_normalized_means, fold_changes_melted, by=c('Representative'))
  } else {
    print("Please provide a valid xlabel")
  }
  
  thresh <- round(quantile(df$value.x, quantile_prob))
  color_vector <- color_vector[unique(df$variable)]
  
  df$value.y <- log2(df$value.y + min(df$value.y[df$value.y!=0]/2))
  p <- ggplot(df, aes(x = value.x, y = value.y)) + 
        geom_point(aes(x = value.x, y = value.y), size = point_size, alpha = 0.7, col="gray35") + 
        geom_point(data=df[((df$value.x >= thresh) & (df$value.y > 1 | df$value.y < -1)),],
                   aes(x = value.x, y = value.y, colour = variable),size=point_size) +
        scale_color_manual(values = color_vector, name="Samples") + 
        geom_hline(yintercept = 1, linetype='dashed', col = 'red', size = 0.4) + 
        geom_hline(yintercept = 0, col = 'red', size = 0.2) + 
        geom_hline(yintercept = -1, linetype='dashed', col = 'red', size = 0.4) + 
        geom_vline(xintercept = thresh, linetype='dashed', col = 'red') +
        coord_cartesian(clip = "off") +
        geom_text_repel(data=df[((df$value.x >= thresh) & (df$value.y > 1 | df$value.y < -1 )),],
                        aes(value.x, value.y,label = Representative, hjust = hjust_), 
                        force = 1, max.overlaps = Inf,
                        min.segment.length = 0.3,
                        xlim = c(-7, Inf), ylim = c(-Inf, Inf),
                        segment.size = 0.1, segment.alpha = 0.6,
                        box.padding = unit(0.25, "lines"),
                        direction = "both",
                        size = 2.2, col = "black") +
        geom_text(aes(thresh,-Inf,label = paste0(thresh, ""), vjust = -0.3, hjust = -0.2), col = "red") +
        ggtitle(title) + 
        xlab(xlabel) + 
        ylab("Log2 fold change") + 
        theme_classic() 
  print(p)
}


