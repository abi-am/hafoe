#Install/Load packages
packages = c("stringr", "seqinr", "microseq")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = T))) {
      install.packages(x, dependencies = TRUE)
      suppressPackageStartupMessages(library(x, character.only = TRUE, quietly = T))
    }
  }
)

# Run after sequencing simulation
for (d in 1:5){
  chimeric_ccs_fq0 <- microseq::readFastq(paste0("input_files/data", d, "/Chimeric_lib_simulated_ccs.fastq.gz"))
  chimeric_ccs_fq <- chimeric_ccs_fq0
  chimeric_ccs_fq$Header <- sub('.*\\=', '', stringr::str_split(chimeric_ccs_fq0$Header, pattern = ";", simplify = TRUE)[,4])
  
  chimeric_enriched_counts <- as.data.frame(table(sub('\\_.*', '', chimeric_ccs_fq$Header)))
  colnames(chimeric_enriched_counts) <- c("name", "chimeric_count")
  
  for (chimseq_name in chimeric_enriched_counts$name){
    chimeric_ccs_fq[chimeric_ccs_fq$Header == chimseq_name,"Name"] <- paste0(chimeric_ccs_fq[chimeric_ccs_fq$Header == chimseq_name,]$Header,
                                                                             "_",
                                                                             1:chimeric_enriched_counts[chimeric_enriched_counts$name == chimseq_name,]$chimeric_count)
  }
  
  chimeric_ccs_fq$Header <- chimeric_ccs_fq$Name
  chimeric_ccs_fq$Name <- NULL
  
  
  chimeric_ccs_csv <- chimeric_ccs_fq[c("Header", "Sequence")]
  colnames(chimeric_ccs_csv) <- c("X", "Sequence")
  chimeric_ccs_csv$Count <- 1
  write.csv(chimeric_ccs_csv, paste0("input_files/data", d, "/Chimeric_lib_simulated_ccs.csv"), row.names = F, quote = F)
}
