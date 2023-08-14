# True vs predicted comporition for synthetic data
##################################################
# PREDICTED

# Multiple sequence alignment of top 20 chimeric library representative
# Output path: file.path(output.dir, "files/top20_rep_msa.clustal_num")

# Inputs
output.dir <- "../../../../hafoe/data_v4/benchmarking/hafoe_out_1_100_10_sum"
simulated_dir <- "../../../../hafoe/data_v4/input_files/data1"
scripts.dir <- "../../src"
clustalo.path <- "/usr/bin/clustalo"
step_size <- 10

source(file.path(scripts.dir, "functions.R"))

cluster_chim_new <- microseq::readFasta(file.path(output.dir, "files/clstr_chimeric_lib/clstr_chim_new.fasta"))
counts_chim <- read.csv(file.path(output.dir, "files/chimeric_lib_rep_counts.csv"))

top20_order <- counts_chim[1:20,]$Representative
top20_order <- top20_order[!is.na(top20_order)]

top20_rep_fa <- cluster_chim_new[gsub("\\/.*","",cluster_chim_new$Header) %in% counts_chim[1:20,]$Representative,]
top20_rep_fa$Header <- gsub("\\/.*","",top20_rep_fa$Header)
top20_rep_fa <- top20_rep_fa[match(top20_order, top20_rep_fa$Header),]

microseq::writeFasta(top20_rep_fa, file.path(output.dir, "files/top20_representatives.fasta"))

o <- system(paste0(file.path(scripts.dir, "msa.sh"),
                   " -c ", clustalo.path,
                   " -i ", file.path(output.dir, "files/top20_representatives.fasta"),
                   " -o ", file.path(output.dir, "files/top20_rep_msa.clustal_num")),
            intern = TRUE) 
cat(o, sep = "\n\n")

# Get conserved and gap regions from multiple sequence alignment
# Output path: file.path(output.dir, "files/variant_description/chimeric_lib_representatives/chimeric_lib_representatives_variant_description_top20_msa.csv"

out <- get.conserved.positions(aln_file_path = file.path(output.dir, "files/top20_rep_msa.clustal_num"))
identity_ranges_top20 <- out[[1]]
alignment_top20 <- out[[2]]


reps_nj_mat <- read.table(file.path(output.dir,
                                    "files/variant_description/chimeric_lib_representatives/chimeric_lib_representatives_variant_description.csv"))
rownames(reps_nj_mat) <- reps_nj_mat[,1]
reps_nj_mat[,1] <- NULL
reps_nj_mat <- as.matrix(reps_nj_mat)
colnames(reps_nj_mat) <- NULL


matrix_nt_top20 <- add.gap.info(alignment = alignment_top20, 
                                nj_matrix = reps_nj_mat, 
                                step_size = step_size, 
                                nj_by_nt = T) 

write.table(matrix_nt_top20, 
            file.path(output.dir, 
                      "files/variant_description", 
                      "chimeric_lib_representatives", 
                      paste0("chimeric_lib_representatives", "_variant_description_top20_msa.csv")), 
            quote = F, 
            col.names = F)

#TRUE
predicted_labels <- read.csv(file.path(output.dir, "files/Chimeric_rep_predicted_labels.csv"))
chim_simulated_initial_fa <- microseq::readFasta(file.path(simulated_dir, "Chimeric_lib_simulated_initial.fasta"))
top20_rep_fa_initial <- chim_simulated_initial_fa[chim_simulated_initial_fa$Header %in% gsub("_.*","",top20_rep_fa$Header),]
top20_rep_fa_initial <-  top20_rep_fa_initial[match(gsub("_.*","",top20_order), top20_rep_fa_initial$Header),]

sim_nj_mat <- read.csv(file.path(simulated_dir,
                                 "Chimeric_lib_simulated_labels.csv"))

sim_nj_mat <- sim_nj_mat[sim_nj_mat$X %in% top20_rep_fa_initial$Header,]
sim_nj_mat[,3] <- NULL
rownames(sim_nj_mat) <- sim_nj_mat[,1]
sim_nj_mat[,1] <- NULL
sim_nj_mat <- as.matrix(sim_nj_mat)
colnames(sim_nj_mat) <- NULL

s <- stringr::str_split(unlist(sim_nj_mat[,1], 1), " ")
col_num <- max(nchar(top20_rep_fa_initial$Sequence))
for (i in seq_len(length(s))){
  if (length(s[[i]]) < col_num){ 
    s[[i]] <- c(s[[i]], rep("18", col_num - length(s[[i]]))) #gap
  }
}

m_new <- matrix(as.numeric(unlist(s)), ncol = col_num, byrow = TRUE)
rownames(m_new) <- rownames(sim_nj_mat)


for (name_i in top20_rep_fa$Header){
  pred <- as.numeric(unlist(stringr::str_split(predicted_labels[predicted_labels$X == name_i,]$Composition, " ")))
  pred <- pred[pred != 18]
  pred_start_orf <- predicted_labels[predicted_labels$X == name_i,]$Start_orf
  pred_end_orf <- predicted_labels[predicted_labels$X == name_i,]$End_orf
  m_new[sub("\\_.*", "", name_i),][1:length(pred)] <- m_new[sub("\\_.*", "", name_i),][pred_start_orf:pred_end_orf]
  m_new[sub("\\_.*", "", name_i),][(length(pred)+1):length(m_new[sub("\\_.*", "", name_i),])] <- 18
  
  seq <- top20_rep_fa_initial[top20_rep_fa_initial$Header == sub("\\_.*", "", name_i), "Sequence"]
  top20_rep_fa_initial[top20_rep_fa_initial$Header == sub("\\_.*", "", name_i), "Sequence"] <- substring(seq,pred_start_orf,pred_end_orf)
}


microseq::writeFasta(top20_rep_fa_initial, file.path(simulated_dir, "Chimeric_lib_simulated_initial_orf_top20.fasta"))

o <- system(paste0(file.path(scripts.dir, "msa.sh"),
                   " -c ", clustalo.path,
                   " -i ", file.path(simulated_dir, "Chimeric_lib_simulated_initial_orf_top20.fasta"),
                   " -o ", file.path(simulated_dir, "Chimeric_lib_simulated_initial_orf_top20.clustal_num")),
            intern = TRUE) 
cat(o, sep = "\n\n")

# Get conserved and gap regions from multiple sequence alignment

out <- get.conserved.positions(aln_file_path = file.path(simulated_dir, "Chimeric_lib_simulated_initial_orf_top20.clustal_num"))
identity_ranges_top20 <- out[[1]]
alignment_top20 <- out[[2]]



matrix_nt_top20_sim <- add.gap.info(alignment = alignment_top20, 
                                    nj_matrix = m_new, 
                                    step_size = step_size, 
                                    nj_by_nt = T) 

write.table(matrix_nt_top20_sim, 
            file.path(simulated_dir,
                      "Chimeric_lib_simulated_labels_top20_msa.csv"), 
            quote = F, 
            col.names = F)

