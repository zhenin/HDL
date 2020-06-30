library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
fn <- gsub(x = args[grep(x = args, pattern = "gwas.file=")], pattern = "gwas.file=", replacement = "")
LD.path <- gsub(x = args[grep(x = args, pattern = "LD.path=")], pattern = "LD.path=", replacement = "")
GWAS.type <- gsub(x = args[grep(x = args, pattern = "GWAS.type=")], pattern = "GWAS.type=", replacement = "")
output.file <- gsub(x = args[grep(x = args, pattern = "output.file=")], pattern = "output.file=", replacement = "")
log.file <- gsub(x = args[grep(x = args, pattern = "log.file=")], pattern = "log.file=", replacement = "")

SNP <- gsub(x = args[grep(x = args, pattern = "SNP=")], pattern = "SNP=", replacement = "")
A1 <- gsub(x = args[grep(x = args, pattern = "A1=")], pattern = "A1=", replacement = "")
A2 <- gsub(x = args[grep(x = args, pattern = "A2=")], pattern = "A2=", replacement = "")
N <- gsub(x = args[grep(x = args, pattern = "N=")], pattern = "N=", replacement = "")
b <- gsub(x = args[grep(x = args, pattern = "b=")], pattern = "b=", replacement = "")
se <- gsub(x = args[grep(x = args, pattern = "se=")], pattern = "se=", replacement = "")
Z <- gsub(x = args[grep(x = args, pattern = "Z=")], pattern = "Z=", replacement = "")

if(length(output.file) == 0)
  output.file <- fn

if(length(log.file) != 0){
  log.file <- paste(log.file, "txt", sep = ".")
  if(file.exists(log.file) == T){
    system(paste0("rm ",log.file))
  }
}




library(data.table)
data.table.version <- packageVersion("data.table")
if(data.table.version < "1.12.1"){
  message("Searching for the updated version of package 'data.table' in user's library:")
  detach("package:data.table", unload=TRUE)
  library(data.table, lib.loc = Sys.getenv("R_LIBS_USER"))
}


smart.reader <- function(path){
  path.split <- unlist(strsplit(path, split = "\\."))
  file.type <- path.split[length(path.split)]
  if(file.type == "rds"){
    return(readRDS(path))
  } else if(file.type == "gz" | file.type == "bgz"){
    options(datatable.fread.input.cmd.message=FALSE)
    return(fread(input = paste("zcat < ",path)))
  } else{
    try_error <- try(return(fread(path)))
    if(class(try_error) == "try-error"){
      error.message <- "This file type is not supported by fread function in data.table package. Please reformat it to .txt, .csv or .tsv."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
  }
}

time.start <- date()
cat("Program starts on",time.start,"\n")
cat("Loading GWAS summary statistics from",fn,"\n")

if(length(log.file) != 0){
  cat("Program starts on",time.start,"\n", file = log.file, append = T)
  cat("Loading GWAS summary statistics from",fn,"\n", file = log.file, append = T)
}
gwas.all <- smart.reader(fn)

cat("Data are loaded successfully. Data wrangling starts. \n")
if(length(log.file) != 0){
  cat("Data are loaded successfully. Data wrangling starts. \n", file = log.file, append = T)
}

LD.files <- list.files(LD.path)
if(any(grepl(x = LD.files, pattern = "UKB_snp_counter.*"))){
  snp_counter_file <- LD.files[grep(x = LD.files, pattern = "UKB_snp_counter.*")]
  snp_list_file <- LD.files[grep(x = LD.files, pattern = "UKB_snp_list.*")]
  load(file=paste(LD.path, snp_counter_file, sep = "/"))
  load(file=paste(LD.path, snp_list_file, sep = "/"))
  if("nsnps.list.imputed" %in% ls()){
    snps.name.list <- snps.list.imputed.vector
    nsnps.list <- nsnps.list.imputed
  }
} else{
  error.message <- "It seems this directory does not contain all files needed for HDL. Please check your LD.path again. The current version of HDL only support pre-computed LD reference panels."
  if(output.file != ""){
    cat(error.message, file = output.file, append = T)
  }
  stop(error.message)
}

## the Neale's UKB GWAS format ##
if(length(GWAS.type) != 0){
  if(GWAS.type == "UKB.Neale"){
    dictionary_file <- LD.files[grep(x = LD.files, pattern = "snp.dictionary.*")]
    load(file=paste(LD.path, dictionary_file, sep = "/"))
    gwas.hdl.df <- gwas.all %>%
      inner_join(snp.dictionary %>% filter(rsid %in% snps.name.list), by = "variant")  %>%
      select(rsid, alt, ref, n_complete_samples, tstat) %>%
      rename(SNP = rsid, A1 = alt, A2 = ref, N = n_complete_samples, Z = tstat)
  }
}

## non built-in format
if(length(GWAS.type) == 0){
  if(length(Z) == 0 & (length(b) == 0 | length(se) == 0)){
    error.message <- "Z-score is not available and either b or se is missing. Please check."
    if(length(log.file) != 0){
      cat(error.message, file = log.file, append = T)
    }
    stop(error.message)
  }

  if(length(Z) != 0){
    gwas.hdl.df <- gwas.all %>%
      rename_(SNP = SNP, A1 = A1, A2 = A2, N = N, Z = Z) %>%
      filter(SNP %in% snps.name.list)
  } else{
    gwas.hdl.df <- gwas.all %>%
      rename_(SNP = SNP, A1 = A1, A2 = A2, N = N, b = b, se = se) %>%
      filter(SNP %in% snps.name.list) %>% mutate(Z = (b/se)) %>%
      select(SNP, A1, A2, N, Z)
  }
}
cat("Data wrangling completed. \n")
if(length(log.file) != 0){
  cat("Data wrangling completed. \n", file = log.file, append = T)
}

k1 <- sum(gwas.hdl.df$SNP %in% snps.name.list)
k1.percent <- paste("(",round(100*k1 / length(snps.name.list), 2), "%)", sep="") 
cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in GWAS."," \n")
if(length(log.file) != 0){
  cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in GWAS."," \n", file = log.file, append = T)
}
if(k1 < length(snps.name.list)*0.99){
  cat("More than 1% SNPs in reference panel are missed in GWAS. This missing rate is too high to be used in HDL. Please check. \n")
  if(length(log.file) != 0){
    cat("More than 1% SNPs in reference panel are missed in GWAS. This missing rate is too high to be used in HDL. Please check. \n", file = log.file, append = T)
  }
}

gwas.hdl.df$A1 <- toupper(gwas.hdl.df$A1)
gwas.hdl.df$A2 <- toupper(gwas.hdl.df$A2)

fn.rds <- paste0(output.file, ".hdl.rds")
saveRDS(gwas.hdl.df, fn.rds)
cat("The output is saved to", fn.rds, "\n")
if(length(log.file) != 0){
  cat("The output is saved to", fn.rds, "\n", file = log.file, append = T)
}

if(length(log.file) != 0){
  cat("The log is saved to", log.file, "\n")
  cat("The log is saved to", log.file, "\n", file = log.file, append = T)
}

