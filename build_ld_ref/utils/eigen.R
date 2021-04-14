#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly=FALSE)
script.name <- sub('--file=', '', initial.options[grep('--file=', initial.options)])
utils.dir <- dirname(script.name)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1){
  # message('Rscript ', script.name, ' <split_window.ld.gz>')
  message('Rscript ', script.name, ' <split_window.ld>')
  stop('command args error.')
}
ld.file <- args[1]



suppressWarnings(suppressMessages(require(RSpectra)))
suppressWarnings(suppressMessages(require(tidyr)))
suppressWarnings(suppressMessages(require(dplyr)))
suppressWarnings(suppressMessages(require(data.table)))

dyn.load(paste0(utils.dir, '/ldscore.so'))
dyn.load(paste0(utils.dir, '/bmult.so'))
source(paste0(utils.dir, '/Rbmult.r'))

# bim.file <- gsub('.ld.gz$', '.bim', ld.file)
# eigen.res.file <- gsub('.ld.gz$', '.banded.rda', ld.file)
bim.file <- gsub('.ld$', '.bim', ld.file)
eigen.res.file <- gsub('.ld$', '.banded.rda', ld.file)
snps <- pull(fread(bim.file, select=2))
bandedR <- fread(ld.file)
nb <- nrow(filter(bandedR, SNP_A==bandedR$SNP_A[1]))
nb <- bandedR %>%
  count(SNP_A, name='count', sort=T) %>%
  head(1) %>%
  select(count) %>%
  pull()

M <- length(snps)
pout <- pull(select(bandedR, R))
nval <- round(length(snps)*0.5)
eigR <- eigs_sym(Rxfun, nval, n=M, which='LA', args=list(nb=nb, pout=pout))
lam <- eigR$values
V <- eigR$vectors
rownames(V) <- snps
LDsc <- Rldscore(M, nb, pout)
save(LDsc, lam, V, file=eigen.res.file)
