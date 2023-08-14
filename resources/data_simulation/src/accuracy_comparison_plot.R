# Run this after running hafoe benchmarking to get accuracy comparison plot

# First specify the packages of interest
packages = c("stringr", "ggplot2")


package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = T))) {
      install.packages(x, dependencies = TRUE)
      suppressPackageStartupMessages(library(x, character.only = TRUE, quietly = T))
    }
  }
)


library(ggplot2)
library(stringr)

# Takes predicted and true labels and calculates the accuracy
get.accuracy <- function(predicted_labels, chimeric_true_labels){
  accuracies <- c()
  for (i in predicted_labels$X) {
    # print(i)
    pred <- as.numeric(unlist(stringr::str_split(predicted_labels[predicted_labels$X == i,]$Composition, " ")))
    pred <- pred[pred != 18]
    pred_start_orf <- predicted_labels[predicted_labels$X == i,]$Start_orf
    pred_end_orf <- predicted_labels[predicted_labels$X == i,]$End_orf
    
    true <- as.numeric(unlist(stringr::str_split(chimeric_true_labels[chimeric_true_labels$X == sub("\\_.*", "", i),]$Composition, " ")))
    #true has labels for the whole sequence, we need only for orf
    true <- true[pred_start_orf:pred_end_orf]  
    
    acc <- round(sum(pred == true)*100/length(true), 2)
    accuracies <- c(accuracies, acc)
  }
  avg_accuracy <- round(sum(accuracies)/length(accuracies), 3)
  return(avg_accuracy)
}


#analyse errors 
get.error.stats <- function(predicted_labels, chimeric_true_labels){
  no_alns <- c()
  multi_alns <- c()
  wrong_alns <- c()
  for (i in predicted_labels$X) {
    pred <- as.numeric(unlist(stringr::str_split(predicted_labels[predicted_labels$X == i,]$Composition, " ")))
    pred <- pred[pred != 18]
    pred_start_orf <- predicted_labels[predicted_labels$X == i,]$Start_orf
    pred_end_orf <- predicted_labels[predicted_labels$X == i,]$End_orf
    
    true <- as.numeric(unlist(stringr::str_split(chimeric_true_labels[chimeric_true_labels$X == sub("\\_.*", "", i),]$Composition, " ")))
    #true has labels for the whole sequence, we need only for orf
    true <- true[pred_start_orf:pred_end_orf]  
    
    no_aln <- sum(pred == 0)*100/length(true)
    multi_aln <- sum(pred == 17)*100/length(true)
    wrong_aln <- sum(((pred != true)&(pred != 0)&(pred != 17)))*100/length(true)
    
    no_alns <- c(no_alns, no_aln)
    multi_alns <- c(multi_alns, multi_aln)
    wrong_alns <- c(wrong_alns, wrong_aln)
  }
  
  print(paste0("No alignments: ", round(mean(no_alns),3), "%"))
  print(paste0("Multiple alignments: ", round(mean(multi_alns), 3), "%"))
  print(paste0("Wrong alignments: ", round(mean(wrong_alns), 3), "%"))

  out <- list(round(mean(no_alns),3), round(mean(multi_alns),3), round(mean(wrong_alns), 3))
  return(out)
}


#get stats for all cases
# count <- 0
columns <- c("data", "rl", "ss", "criterion", "accuracy", "no_alignment", "multiple_alignment", "wrong_alignment") 
acc_df <- data.frame(matrix(nrow=1, ncol = length(columns))) 
colnames(acc_df) <- columns

for (d in 1:5){
  for (rl in seq(50, 200, 50)){
    for (ss in seq(10, 200, 10)){
      if (rl >= ss){
        for (criterion in c("avg", "sum")){
          path_pred <- paste0("benchmarking/hafoe_out_", d, "_", rl, "_", ss, "_", criterion, "/files/Chimeric_rep_predicted_labels.csv")
          path_true <- paste0("input_files/data", d, "/Chimeric_lib_simulated_labels.csv")
          predicted_labels_temp <- read.csv(path_pred)
          true_labels_temp <- read.csv(path_true)
          acc <- get.accuracy(predicted_labels_temp, true_labels_temp)
          out <- get.error.stats(predicted_labels_temp, true_labels_temp)
          no_aln <- out[[1]]
          mult_aln <- out[[2]]
          wrong_aln <- out[[3]]
          
          acc_df <- rbind(acc_df, c(d, rl, ss, criterion, acc, no_aln, mult_aln, wrong_aln))
          # if(is.na(a)){
          #   print(path_pred) 
          #   count <- count + 1
          # }
        }
      }
    }
  }
}

acc_df <- na.omit(acc_df)
acc_df$accuracy <- as.numeric(acc_df$accuracy)
acc_df$rl <- as.numeric(acc_df$rl)
acc_df$ss <- as.numeric(acc_df$ss)
acc_df <- acc_df[order(as.numeric(acc_df$ss), decreasing = T),]
acc_df <- acc_df[order(as.numeric(acc_df$rl)),]


# PLot accuracy comaprison boxplot

# facet label names 
rl.labs <- paste0(unique(acc_df$rl), "nt read length")
names(rl.labs) <- unique(acc_df$rl)

pdf(file.path("plots/variant_composition_accuracy_comparison.pdf"), width=11, height=6)
ggplot(acc_df, aes(x=factor(ss, levels = sort(unique(as.numeric(ss)), decreasing = T)), 
                   y=accuracy, 
                   fill=criterion)) + 
  geom_boxplot(lwd=0.4) +
  scale_x_discrete(drop = T) + 
  facet_wrap(~rl, scales = "free_x", labeller = labeller(rl = rl.labs)) + 
  labs(title=paste0(""),
       x ="Step size, nt", y = "Accuracy (%)") +
  theme_bw()
dev.off()