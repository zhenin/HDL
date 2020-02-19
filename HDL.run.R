args <- commandArgs(trailingOnly = TRUE)
gwas1.df.path <- gsub(x = args[grep(x = args, pattern = "gwas1.df=")], pattern = "gwas1.df=", replacement = "")
gwas2.df.path <- gsub(x = args[grep(x = args, pattern = "gwas2.df=")], pattern = "gwas2.df=", replacement = "")
LD.path <- gsub(x = args[grep(x = args, pattern = "LD.path=")], pattern = "LD.path=", replacement = "")
Nref <- gsub(x = args[grep(x = args, pattern = "Nref=")], pattern = "Nref=", replacement = "")
N0 <- gsub(x = args[grep(x = args, pattern = "N0=")], pattern = "N0=", replacement = "")
output.file <- gsub(x = args[grep(x = args, pattern = "output.file=")], pattern = "output.file=", replacement = "")


# gwas1.df.path <- "/Users/zhengning/Work/HDL/package/HDL_examples/gwas1.example.txt"
# gwas2.df.path <- "/Users/zhengning/Work/HDL/package/HDL_examples/gwas2.example.rds"
smart.reader <- function(path){
  path.split <- unlist(strsplit(path, split = "\\."))
  file.type <- path.split[length(path.split)]
  if(file.type == "rds"){
    return(readRDS(path))
  } else if(file.type == "txt"){
    return(read.table(file = path, header = T))
  } else{
    stop("The extension of input file has to be .rds or .txt!")
  }
}
gwas1.df <- smart.reader(gwas1.df.path)
gwas2.df <- smart.reader(gwas2.df.path)

if(length(Nref)==0)
  Nref <- 336000
if(length(N0) == 0)
  N0 <- min(gwas1.df$N, gwas2.df$N)
if(length(output.file) == 0)
  output.file <- ""



##### Run HDL #####
if(output.file != ""){
  if(file.exists(output.file) == T){
    file.remove(output.file)
  }
  capture.output(library(HDL), type = "message", file = output.file)
  capture.output(res.HDL <- HDL.rg(gwas1.df, gwas2.df, LD.path, Nref = Nref, N0 = N0, output.file = output.file), 
                 file = output.file, append = TRUE)
} else{
  library(HDL)
  res.HDL <- HDL.rg(gwas1.df, gwas2.df, LD.path, Nref = Nref, N0 = N0, output.file = output.file)
}

if(abs(res.HDL$rg) < 1e-4){
  rg.out <- formatC(res.HDL$rg, format = "e", digits = 2)
} else{
  rg.out <- round(res.HDL$rg, digits = 4)
}

if(res.HDL$rg.se < 1e-4){
  rg.se.out <- formatC(res.HDL$rg.se, format = "e", digits = 2)
} else{
  rg.se.out <- round(res.HDL$rg.se, digits = 4)
}

p.out <- formatC(res.HDL$P, format = "e", digits = 2)

cat("\n", file = output.file, append = TRUE)
cat("Genetic Correlation: ", 
    rg.out, 
    paste0("(", rg.se.out, ")"), 
    file = output.file, append = TRUE)
cat("\n", file = output.file, append = TRUE)
cat("P: ",p.out, file = output.file, append = TRUE)
cat("\n", file = output.file, append = TRUE)
cat("\n", file = output.file, append = TRUE)
cat("Analysis finished at",date(),"\n", file = output.file, append = TRUE)
if(output.file != ""){
  cat("The results were saved to", output.file, file = output.file, append = TRUE)
}
cat("\n", file = output.file, append = TRUE)