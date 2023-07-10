args <- commandArgs(TRUE)

config.file = args[1]


if (length(args) < 1){  
  stop("No config file provided")
}

if (file.exists(config.file))
  config.set = T else stop("No valid config file provided")

#Read the options from config file
argnames = vector()
argvalues = vector()
lines = readLines(config.file)
for(l in 1:length(lines)){
  line = lines[l]
  line = unlist(strsplit(line, split = "#"))[1]
  line = unlist(strsplit(line, split = "\\*"))[1]  
  tokens = unlist(strsplit(line, split = " ")) #\t
  if(length(tokens) == 2){
    argvalues = c(argvalues, tokens[2])
    argnames = c(argnames, tokens[1])
  }
}
config.table = matrix(argvalues, length(argvalues), 1)
rownames(config.table) = argnames


# config.table = as.matrix(read.table(file=config.file, header=F, sep="\t", comment.char =c("#","*"), blank.lines.skip=T, as.is=TRUE, row.names=1))

#######################################################################
##########	Validate the options	###############
#######################################################################


set.property <- function(matrix, prop.name, default.value){
  result = tryCatch({
    matrix[prop.name,1]
  }, warning = function(w) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  })
  return(result)
}

set.property.double <- function(matrix, prop.name, default.value){
  result = set.property(matrix, prop.name, default.value)
  result = tryCatch({
    as.double(result)
  }, warning = function(w) {
    cat("Warning: Illegal value for", prop.name, ": ", result,
        "Should be double. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("Error: Illegal value for", prop.name, ": ", result,
        "Should be double. Default value assigned: ", default.value, "\n")
    default.value
  })
}

set.property.integer <- function(matrix, prop.name, default.value){
  result = set.property(matrix, prop.name, default.value)
  result = tryCatch({
    as.integer(result)
  }, warning = function(w) {
    cat("Warning: Illegal value for", prop.name, ": ", result,
        "Should be integer. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("Error: Illegal value for", prop.name, ": ", result,
        "Should be integer. Default value assigned: ", default.value, "\n")
    default.value
  })
}
set.property.logical <- function(matrix, prop.name, default.value){
  result = set.property(matrix, prop.name, default.value)
  result = tryCatch({
    as.logical(result)
  }, warning = function(w) {
    cat("Warning: Illegal value for", prop.name, ": ", result,
        "Should be logical. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("Error: Illegal value for", prop.name, ": ", result,
        "Should be logical. Default value assigned: ", default.value, "\n")
    default.value
  })
}


set.property.executable <- function(matrix, prop.name, default.value){
  result = tryCatch({
    matrix[prop.name,1]
  }, warning = function(w) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  }, error= function(e) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  })
  if (length(grep("./", result, fixed = T)) > 0
      || length(grep(".\\", result, fixed = T)) > 0){
    wd = getwd()
    setwd("..")
    parent = getwd()
    setwd(wd)
  }
  
  if (length(grep(pattern="../", result, fixed = T)) > 0){
    result = sub(pattern="..", replacement=parent, result, fixed = T)
  } else if (length(grep(pattern="./", result, fixed = T)) > 0){
    result = sub(pattern=".", replacement=wd, result, fixed = T)
  } else   if (length(grep(pattern="..\\", result, fixed = T)) > 0){
    result = sub(pattern="..", replacement=parent, result, fixed = T)
  } else if (length(grep(pattern=".\\", result, fixed = T)) > 0){
    result = sub(pattern=".", replacement=wd, result, fixed = T)
  }
  
  return(result)
}

defaults = list(
  output.dir = "hafoe_out",
  explore = F,
  identify = F,
  read_length = 100,
  step_size = 10,
  vd_criterion = "sum",
  scripts.dir = "./",
  bowtie.dir = "../bowtie2-2.1.0/",
  bowtie.build.path = "../bowtie2-2.1.0/bowtie2-build.exe",
  bowtie.align.path = "../bowtie2-2.1.0/bowtie2-align.exe",
  samtools.path = "../samtools-0.1.19/samtools.exe",
  ignore.err = F
)


explore = set.property.logical(config.table, "explore", defaults$explore)
identify = set.property.logical(config.table, "identify", defaults$identify)
if(!explore & !identify){
  cat("Error: explore and identify cannot be false at the same time\n")
  config.set = F
}

read_length = set.property.integer(config.table, "read_length", defaults$read_length)
step_size = set.property.integer(config.table, "step_size", defaults$step_size)

vd_criterion = set.property(config.table, "vd_criterion", defaults$vd_criterion)



parent.path = tryCatch({
  config.table["parent",1]
}, warning = function(w) {
  cat("\"parent\" parameter not specified\n")
  config.set = F
}, error = function(e) {
  cat("\"parent\" parameter not specified\n")
  config.set = F
})

chimeric.lib.path = tryCatch({
  config.table["chimera",1]
}, warning = function(w) {
  cat("\"chimera\" parameter not specified\n")
  config.set = F
}, error = function(e) {
  cat("\"chimera\" parameter not specified\n")
  config.set = F
})

explore.out.path = tryCatch({
  config.table["explore.out",1]
}, warning = function(w) {
  cat("\"explore.out\" parameter not specified\n")
  config.set = F
}, error = function(e) {
  cat("\"explore.out\" parameter not specified\n")
  config.set = F
})

enriched.path = tryCatch({
  config.table["enriched",1]
}, warning = function(w) {
  cat("\"enriched\" parameter not specified\n")
  config.set = F
}, error = function(e) {
  cat("\"enriched\" parameter not specified\n")
  config.set = F
})
names(enriched.path) <- NULL

# if (enriched == "none"){
#   cat("no enrichment file provided.")
# } else {
if (enriched.path != "none"){
  enriched.path = unlist(strsplit(enriched.path, split='[ ,]' , fixed=F))
  enriched.path = enriched.path[!duplicated(enriched.path)]
  empty.enriched = c()
  for (i in 1:length(enriched.path)){
    if (identical(enriched.path[i], "")) {
      empty.enriched = c(empty.enriched,i)
    } else if (!file.exists(enriched.path[i])){
      cat("en file ", enriched.path[i],
          " does not exist. Please, specify a valid file.\n")
      config.set = F
    }
  }
  if(length(empty.enriched) > 0)
    enriched.path=enriched.path[-empty.enriched]
  
  en1.path = tryCatch({
    enriched.path[1]
  }, warning = function(w) {
    cat("\"enriched1\" parameter not specified\n")
    config.set = F
  }, error = function(e) {
    cat("\"enriched1\" parameter not specified\n")
    config.set = F
  })
  
  #input multiple files
  en1.paths <- c()
  for (en1 in sort(list.files(en1.path, full.names=T))){
    if (grepl("fastq", en1) | grepl("fq", en1)){
      en1.paths <- c(en1.paths, en1)
    } else {
      stop(paste0("no enriched1 fastq files found at ", en1.path, ".\n"))
    }
  }

  if (length(en1.paths) == 0){
    stop(paste0("no enriched1 files found at ", en1.path, ".\n"))
  }

  
  
  en2.path = tryCatch({
    enriched.path[2]
  }, warning = function(w) {
    cat("\"enriched2\" parameter not specified\n")
    config.set = F
  }, error = function(e) {
    cat("\"enriched2\" parameter not specified\n")
    config.set = F
  })
  
  #input multiple files
  # en2.paths <- c()
  # for (en2 in sort(list.files(en2.path, full.names=T))){
  #   if (grepl("fastq", en2) | grepl("fq", en1)){
  #     en2.paths <- c(en2.paths, en2)
  #   } else {
  #     stop(paste0("no enriched2 fastq files found at ", en2.path, ".\n"))
  #   }
  # }
  # 
  # if (length(en2.paths) == 0){
  #   stop(paste0("no enriched2 files found at ", en2.path, ".\n"))
  # }
  # 
  # en1.path <- enriched.path[1]
  # en2.path <- enriched.path[2]
}




output.dir = set.property(config.table, "output.dir", defaults$output.dir)
if(is.na(file.info(output.dir)$isdir || !file.info(output.dir)$isdir)){
  dc = dir.create(output.dir, showWarnings=F)
  if (!dc){
    cat("Warning: could not create output directory", output.dir,
        "results will be kept in default directory", defaults$output.dir, "\n")
    output.dir = defaults$output.dir
    dir.create(output.dir, showWarnings=F)
  }
}


rlib = set.property(config.table, "rlib", file.path(output.dir, "rlib"))
if(is.na(file.info(rlib)$isdir || !file.info(rlib)$isdir)){
  dc = dir.create(rlib, showWarnings=F)
  if (!dc){
    cat("Warning: could not create output directory", rlib,
        "results will be kept in default directory",  file.path(output.dir, "rlib"), "\n")
    rlib = file.path(output.dir, "rlib")
    dir.create(rlib, showWarnings=F)
  }
}



scripts.dir = set.property(config.table, "scripts.dir", defaults$scripts.dir)
scripts.dir.set = F
if (file.exists(file.path(scripts.dir, "pipeline.R")))
  if (file.exists(file.path(scripts.dir, "functions.R"))){
    scripts.dir.set = T
  }
#}
if (!scripts.dir.set){
  cat("scripts.dir does not exist or does not contain the required scripts\n")
  config.set = F
}


bowtie.build.path = set.property.executable(config.table, "bowtie.build.path", defaults$bowtie.build.path)
if (!file.exists(bowtie.build.path)){
  cat("bowtie.build executable could not be found at ", bowtie.build.path, "\n")
  config.set = F
}

bowtie.align.path = set.property.executable(config.table, "bowtie.align.path", defaults$bowtie.align.path)
if (!file.exists(bowtie.align.path)){
  cat("bowtie.align executable not be found at ", bowtie.align.path, "\n")
  config.set = F
}

samtools.path = set.property.executable(config.table, "samtools.path", defaults$samtools.path)
if (!file.exists(samtools.path)){
  cat("samtools.path", samtools.path, "does not exist\n")
  config.set = F
}

cd_hit_est.path = set.property.executable(config.table, "cd_hit_est.path", defaults$cd_hit_est.path)
if (!file.exists(cd_hit_est.path)){
  cat("cd_hit_est.path", cd_hit_est.path, "does not exist\n")
  config.set = F
}


cd_hit_est_2d.path = set.property.executable(config.table, "cd_hit_est_2d.path", defaults$cd_hit_est_2d.path)
if (!file.exists(cd_hit_est_2d.path)){
  cat("cd_hit_est_2d.path", cd_hit_est_2d.path, "does not exist\n")
  config.set = F
}

clustalo.path = set.property.executable(config.table, "clustalo.path", defaults$clustalo.path)
if (!file.exists(clustalo.path)){
  cat("clustalo.path", clustalo.path, "does not exist\n")
  config.set = F
}


# Create log file
#######################################################################################
#######################################################################################

dir.create(file.path(output.dir, "log"), showWarnings = F)
log_file <- file(file.path(output.dir, "log/hafoe.log"), open = "wt")
sink(log_file, append = T, type = "output") 
sink(log_file, append = T, type = "message") 

if (config.set) {
  
  cat("\nConfiguration successful\n")
  cat("===================================================================\n\n")
  cat("hafoe will execute with the following parameters:\n")
  cat("\texplore =", explore, "\n")
  cat("\tidentify =", identify, "\n")
  cat("\tscripts.dir =", scripts.dir, "\n")
  cat("\tbowtie.build.path =", bowtie.build.path, "\n")
  cat("\tbowtie.align.path =", bowtie.align.path, "\n")
  cat("\tsamtools.path =", samtools.path, "\n")
  cat("\tcd_hit_est.path =", cd_hit_est.path, "\n")
  cat("\tcd_hit_est_2d.path =", cd_hit_est_2d.path, "\n")
  cat("\tclustalo.path =", clustalo.path, "\n")
  cat("\tparent.path =", parent.path, "\n")
  cat("\tchimeric.lib.path =", chimeric.lib.path, "\n")
  cat("\tenriched.path =", enriched.path, "\n")
  cat("\toutput.dir = ", output.dir,"\n")
  cat("\trlib = ", rlib,"\n")
  cat("\tread_length =", read_length, "\n")
  cat("\tstep_size =", step_size, "\n")
  cat("\tvd_criterion =", vd_criterion, "\n")

}

sink()


if (!config.set){
  stop("configuration not set successfully. Scripts will not execute.\n")
} else {
  source(file.path(scripts.dir, "pipeline.R"))
}