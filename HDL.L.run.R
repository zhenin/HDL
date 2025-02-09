library(parallel)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
  library(HDL)
  q()
}
args.print <- paste("Function arguments:", paste(args, collapse = "\n"), sep = "\n")
cat(args.print, "\n\n")
gwas1.df.path <- gsub(x = args[grep(x = args, pattern = "gwas1.df=")], pattern = "gwas1.df=", replacement = "")
gwas2.df.path <- gsub(x = args[grep(x = args, pattern = "gwas2.df=")], pattern = "gwas2.df=", replacement = "")
Trait1name <- gsub(x = args[grep(x = args, pattern = "Trait1name=")], pattern = "Trait1name=", replacement = "")
Trait2name <- gsub(x = args[grep(x = args, pattern = "Trait2name=")], pattern = "Trait2name=", replacement = "")
LD.path <- gsub(x = args[grep(x = args, pattern = "LD.path=")], pattern = "LD.path=", replacement = "")
bim.path <- gsub(x = args[grep(x = args, pattern = "bim.path=")], pattern = "bim.path=", replacement = "")
Nref <- gsub(x = args[grep(x = args, pattern = "Nref=")], pattern = "Nref=", replacement = "")
N0 <- gsub(x = args[grep(x = args, pattern = "N0=")], pattern = "N0=", replacement = "")
output.file <- gsub(x = args[grep(x = args, pattern = "output.file=")], pattern = "output.file=", replacement = "")
eigen.cut <- gsub(x = args[grep(x = args, pattern = "eigen.cut=")], pattern = "eigen.cut=", replacement = "")
intercept.output <- gsub(x = args[grep(x = args, pattern = "intercept.output=")], pattern = "intercept.output=", replacement = "")
fill.missing.N <- gsub(x = args[grep(x = args, pattern = "fill.missing.N=")], pattern = "fill.missing.N=", replacement = "")
lim  <- gsub(x = args[grep(x = args, pattern = "lim=")], pattern = "lim=", replacement = "")
type <- gsub(x = args[grep(x = args, pattern = "type=")], pattern = "type=", replacement = "")
chr <- gsub(x = args[grep(x = args, pattern = "chr=")], pattern = "chr=", replacement = "")
piece <- gsub(x = args[grep(x = args, pattern = "piece=")], pattern = "piece=", replacement = "")
cores <- gsub(x = args[grep(x = args, pattern = "cores=")], pattern = "cores=", replacement = "")
save.path <- gsub(x = args[grep(x = args, pattern = "save.path=")], pattern = "save.path=", replacement = "")
library(HDL)
## default arguments
Nref <- ifelse(length(Nref)==0, 335272, as.numeric(Nref))
N0 <- ifelse(length(N0)==0, 0, as.numeric(N0))
eigen.cut <- ifelse(length(eigen.cut)==0, 0.99, as.numeric(eigen.cut))
if(!any(fill.missing.N==c("median", "min", "max"))) fill.missing.N <- NULL
output.file <- ifelse(length(output.file)==0, "", output.file)
cores <- ifelse(length(cores)==0, 1, as.numeric(cores))
type <- ifelse(length(type)==0, "WG", type)
chr <- ifelse(length(chr)==0, "", as.numeric(chr))
piece <- ifelse(length(piece)==0, "", as.numeric(piece))
lim <- ifelse(length(lim)==0, exp(-18), as.numeric(lim))


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
    return(as.data.frame(readRDS(path)))
  } else if(file.type == "gz" | file.type == "bgz"){
    options(datatable.fread.input.cmd.message=FALSE)
    return(as.data.frame(fread(input = paste("zcat < ",path))))
  } else{
    try_error <- try(return(as.data.frame(fread(path))))
    if(!is.null(try_error)){
      error.message <- "This file type is not supported by fread function in data.table package. Please reformat it to .txt, .csv or .tsv."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
  }
}

if(length(gwas1.df.path) != 0 && length(gwas2.df.path) != 0){
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
  
  if(length(N0) == 0)
    N0 <- min(gwas1.df$N, gwas2.df$N)
  
}

load(paste0(LD.path, "HDLL_LOC_snps.RData"))


process_function <- function(chr, piece) {
  cat("Processing chromosome", chr, "region", piece, "\n")
  tryCatch({
    res <- HDL.L(gwas1.df, gwas2.df, Trait1name = Trait1name, Trait2name = Trait2name,
                 LD.path = LD.path, bim.path = bim.path, chr = chr, piece = piece,
                 output.file = output.file, N0 = N0, lim = lim, Nref = Nref, eigen.cut = eigen.cut)
    return(res)
  }, error = function(e) {
    cat("Error in region", piece, ":", conditionMessage(e), "\n")
    return(data.frame())  # Return an empty data.frame on error to keep the output consistent
  })
}

if(type == "WG"){
    results <- mcmapply(process_function, as.numeric(NEWLOC$CHR), as.numeric(NEWLOC$piece), SIMPLIFY = FALSE, mc.cores = cores)
}else if(length(chr) > 0 && length(piece) > 0 && type =="region"){
    results <- mcmapply(process_function, chr, piece, SIMPLIFY = FALSE, mc.cores = 1)
}

# Combine all results into a single data frame
res.HDLL <- do.call(rbind, results)
cat("Saving results ... \n")
save(res.HDLL, file = paste0(save.path,"res_HDLL",Trait1name,"_", Trait2name,".RData"))


  if(output.file != ""){
    fConn <- file(output.file)
    Lines <- readLines(fConn)
    writeLines(c(args.print,Lines), con = fConn)
    close(fConn)
  }

cat("Finished!\n")
if(output.file != ""){
  message <- "Finished!\n"
  cat(message, file = output.file, append = T)
}
