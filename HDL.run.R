args <- commandArgs(trailingOnly = TRUE)
gwas1.df.path <- gsub(x = args[grep(x = args, pattern = "gwas1.df=")], pattern = "gwas1.df=", replacement = "")
gwas2.df.path <- gsub(x = args[grep(x = args, pattern = "gwas2.df=")], pattern = "gwas2.df=", replacement = "")
LD.path <- gsub(x = args[grep(x = args, pattern = "LD.path=")], pattern = "LD.path=", replacement = "")
Nref <- gsub(x = args[grep(x = args, pattern = "Nref=")], pattern = "Nref=", replacement = "")
N0 <- gsub(x = args[grep(x = args, pattern = "N0=")], pattern = "N0=", replacement = "")
output.file <- gsub(x = args[grep(x = args, pattern = "output.file=")], pattern = "output.file=", replacement = "")


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gwas1.df <- loadRData(gwas1.df.path)
gwas2.df <- loadRData(gwas2.df.path)

if(length(Nref)==0)
  Nref <- 336000
if(length(N0) == 0)
  N0 <- min(gwas1.df$N, gwas2.df$N)
if(length(output.file) == 0)
  output.file <- ""

library(HDL)

##### Run HDL #####
if(output.file != ""){
  capture.output(res.HDL <- HDL.rg(gwas1.df, gwas2.df, LD.path, Nref = Nref, N0 = N0, output.file = output.file), 
                 file = output.file, append = TRUE)
} else{
  res.HDL <- HDL.rg(gwas1.df, gwas2.df, LD.path, Nref = Nref, N0 = N0, output.file = output.file)
}

cat("\n", file = output.file, append = TRUE)
cat("Genetic Correlation: ", round(res.HDL$rg, digits = 4), paste0("(", round(res.HDL$rg.se, digits = 4), ")"), 
    file = output.file, append = TRUE)
cat("\n", file = output.file, append = TRUE)

