#' Re-build the LD reference panel with given SNPs.
#'
#' @param LD.path Path to the \code{.rda} file where the Eigen decomposition of LD matrix is stored.
#' @param gwas.snps A vector contains SNP IDs from the GWAS summary statistics.
#' @param LD.path.new The new path to the re-built LD reference.
#' @param nthreads Number of threads to use, default \code{nthreads = 1}.
#' @param lam.cut Eigenvalue cutoff for LD matrices, default \code{lam.cut = NULL}, which means the same cutoff in the original LD reference panel. For analyses with a limited number of traits and annotations, a lower cutoff (such as 0.1, or even not using a cutoff at all) is recommended. For large-scale analyses, a higher cutoff (such as 1)  is recommended, to yield fast computation.
#' @param pattern Chromosome and picece pattern of LD files, default is \code{".*chr(\\\\d{1,2})\\\\.(\\\\d{1,2})[_\\\\.].*"}.
#' @return The new path to the re-built LD reference.
#' @export
#' 
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads

HDL.S.rebuild.ref <- function(LD.path, gwas.snps, LD.path.new, nthreads=1,
  lam.cut=NULL, pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"){
  # use RhpcBLASctl to control the number of thread
  blas_set_num_threads(nthreads)
  omp_set_num_threads(nthreads)
  LD.files <- HDL.S:::list.LD.ref.files(LD.path, suffix=".rda", pattern=pattern)
  bim.files <- HDL.S:::list.LD.ref.files(LD.path, suffix=".bim", pattern=pattern)
  LD.path.new <- normalizePath(LD.path.new, mustWork=FALSE)
  if(!dir.exists(LD.path.new)) dir.create(LD.path.new, recursive = TRUE)
  # clean the new directory
  HDL.S:::log.msg("Cleaning the new directory...\n")
  files <- list.files(LD.path.new, full.names=TRUE)
  if(length(files) > 0) file.remove(files)

  chroms <- gsub(pattern, "\\1", basename(LD.files))
  chroms <- unique(as.numeric(chroms))
  gwas.snps <- as.character(gwas.snps)
  gwas.snps <- gwas.snps[!duplicated(gwas.snps)]
  nsnps.list <- rep(list(numeric()), length(chroms))
  names(nsnps.list) <- chroms
  snps.name.list <- c()
  for(i in seq_along(LD.files)){
    t0 <- Sys.time()
    chrom <- as.numeric(gsub(pattern, "\\1", basename(LD.files[i])))
    piece <- as.numeric(gsub(pattern, "\\2", basename(LD.files[i])))
    load(LD.files[i]) # V, lam, LDsc
    out.LD <- file.path(LD.path.new, basename(LD.files[i]))
    out.bim <- file.path(LD.path.new, basename(bim.files[i]))
    if(exists('U') && !exists('V')){
      V <- U
    }
    if(is.null(lam.cut)) lam.cut <- min(lam) - 1e-9
    M <- nrow(V)
    kept.snps <- intersect(rownames(V), gwas.snps)
    R <- V %*% diag(lam) %*% t(V)
    colnames(R) <- rownames(R) <- rownames(V)
    R <- R[kept.snps, kept.snps]
    LDsc <- as.numeric(diag(crossprod(R)))
    eig <- eigen(R, symmetric=TRUE)
    lam <- eig$values
    V <- eig$vectors
    idx <- which(lam > lam.cut)
    lam <- lam[idx]
    V <- V[, idx]
    rownames(V) <- kept.snps
    save(V, lam, LDsc, file=out.LD)

    bim <- read.table(bim.files[i], header=FALSE, stringsAsFactors=FALSE)
    bim <- bim[match(kept.snps, bim$V2), ]
    write.table(bim, out.bim, quote=FALSE, row.names=FALSE, col.names=FALSE)
    nsnps.list[[chrom]] <- c(nsnps.list[[chrom]], length(kept.snps))
    snps.name.list <- c(snps.name.list, kept.snps)
    tt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    HDL.S:::log.msg(sprintf(
      "Updating chr%d_%d in %.2f seconds: %d (out of %d) SNPs are kept.\n",
      chrom, piece, tt, length(kept.snps), M
    ))
  }
  save(nsnps.list, file=file.path(LD.path.new, 'updated_snp_counter_array.RData'))
  save(snps.name.list, file=file.path(LD.path.new, 'updated_snp_list_array.RData'))
  return(LD.path.new)
}