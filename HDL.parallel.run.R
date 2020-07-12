if(!require("doSNOW",character.only = TRUE)){
  cat("Pacakge doSNOW has not been installed. Start installing... \n")
  try_error <- try(install.packages("doSNOW",dependencies=TRUE, lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  if(!is.null(try_error)){
    try_error <- try(install.packages("doSNOW", 
                                      dependencies=TRUE), silent = TRUE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
args.print <- paste("Function arguments:", paste(args, collapse = "\n"), sep = "\n")
cat("\n")
cat(args.print, "\n\n")
if(length(args) == 0){
  library(HDL)
  q()
}
gwas1.df.path <- gsub(x = args[grep(x = args, pattern = "gwas1.df=")], pattern = "gwas1.df=", replacement = "")
gwas2.df.path <- gsub(x = args[grep(x = args, pattern = "gwas2.df=")], pattern = "gwas2.df=", replacement = "")
LD.path <- gsub(x = args[grep(x = args, pattern = "LD.path=")], pattern = "LD.path=", replacement = "")
Nref <- gsub(x = args[grep(x = args, pattern = "Nref=")], pattern = "Nref=", replacement = "")
N0 <- gsub(x = args[grep(x = args, pattern = "N0=")], pattern = "N0=", replacement = "")
output.file <- gsub(x = args[grep(x = args, pattern = "output.file=")], pattern = "output.file=", replacement = "")
numCores <- gsub(x = args[grep(x = args, pattern = "numCores=")], pattern = "numCores=", replacement = "")
eigen.cut <- gsub(x = args[grep(x = args, pattern = "eigen.cut=")], pattern = "eigen.cut=", replacement = "")
jackknife.df <- gsub(x = args[grep(x = args, pattern = "jackknife.df=")], pattern = "jackknife.df=", replacement = "")



if(length(output.file) == 0){
  length(output.file) <- ""
}

if(output.file != ""){
  if(file.exists(output.file) == T){
    system(paste0("rm ",output.file))
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
    if(!is.null(try_error)){
      error.message <- "This file type is not supported by fread function in data.table package. Please reformat it to .txt, .csv or .tsv."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
  }
}

message <- "Loading GWAS1 ... \n"
if(output.file != ""){
  cat(message, file = output.file, append = T)
}
cat(message)
gwas1.df <- smart.reader(gwas1.df.path)

message <- "Loading GWAS2 ... \n"
if(output.file != ""){
  cat(message, file = output.file, append = T)
}
cat(message)
gwas2.df <- smart.reader(gwas2.df.path)

if(length(Nref)==0)
  Nref <- 335265
if(length(N0) == 0)
  N0 <- min(gwas1.df$N, gwas2.df$N)
if(length(eigen.cut)==0)
  eigen.cut <- "automatic"
if(length(jackknife.df)==0)
  jackknife.df <- FALSE


##### Run HDL #####

library(HDL)
res.HDL <- HDL.rg.parallel(gwas1.df, gwas2.df, LD.path, Nref = Nref, N0 = N0, output.file = output.file, 
                           numCores = numCores, eigen.cut = eigen.cut, jackknife.df = jackknife.df)

if(output.file != ""){
  fConn <- file(output.file)
  Lines <- readLines(fConn)
  writeLines(c(args.print,Lines), con = fConn)
  close(fConn)
  if(jackknife.df == TRUE){
    write.table(res.HDL$jackknife.df, file = paste0(output.file,"_jackknife.df",".txt"))
  }
}