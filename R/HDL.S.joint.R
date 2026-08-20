#' Maximize the log likelihood to infer the heritability enrichment folds jointly.
#'
#' @param Dr.path The path to the pre-computed \code{diag(Dr)} data, see also \code{\link{HDL.S::HDL.S.joint.reduct.dim}}.
#' @param gwas.df A data frame including GWAS summary statistics of genetic variants for a trait.
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS.
#' @param LD.path Path to the directory where linkage disequilibrium (LD) information is stored, compatible with HDL LD reference panels.
#' @param output.file Where the log and results should be written. If you do not specify a file, the log will be printed on the console.
#' @param annots A vector of annotation names, default \code{NULL}, which means all annotations in the \code{Dr.path} will be used.
#' @param mc.cores Number of cores to use for parallelization, default \code{mc.cores = 1}.
#' @param stepwise Whether to estimate enrichment fold by estimating heritability, intercept and eta first, default \code{stepwise = FALSE}. If \code{fix.h2} and \code{fix.intercept} are not \code{NULL}, \code{stepwise} will be overridden.
#' @param fill.missing.N If the sample size is missing in the GWAS summary statistics, \code{fill.missing.N} is used to fill the missing values.
#' @param lim Tolerance limitation to ensure the positive-definiteness of covariance matrices, default \code{lim = exp(-18)}.
#' @param lam.cut Eigenvalue cutoff for LD matrices, default \code{lam.cut = NULL}, which means no cutoff. For analyses with a limited number of traits and annotations, a lower cutoff (such as 0.1, or even not using a cutoff at all) is recommended. For large-scale analyses, a higher cutoff (such as 1)  is recommended, to yield fast computation.
#' @param verbose Whether to print the log on the console, default \code{verbose = FALSE}.
#' @param fix.h2 Whether to fix the heritability to \code{fix.h2} or estimate the heritability, default \code{fix.h2 = NULL}, which means estimate the heritability.
#' @param fix.intercept Whether to fix the intercept to \code{fix.intercept} or estimate the intercept, default \code{fix.intercept = NULL}, which means estimate the intercept.
#' @param fix.eta Whether to fix the eta to \code{fix.eta} or estimate the eta, default \code{fix.eta = NULL}, which means estimate the eta. If \code{fix.eta} is not \code{NULL}, the eta will be fixed to \code{fix.eta}.
#' @param maxit Maximum number of iterations, default \code{maxit = 1000}.
#' @param pgtol Tolerance for convergence, default \code{pgtol = 1e-3}.
#' @param start.v Starting values for \code{c(h2, intercept, eta, fold1, fold2, ..., foldn)}, where n is the number of annotations. Default is \code{NULL}, which means \code{c(0.1, 1, 0, 1, ..., 1)}. If \code{fix.h2}, \code{fix.intercept} or \code{fix.eta} is not \code{NULL}, the starting values for \code{h2}, \code{intercept} or \code{eta} will be ignored.
#' @param lwr Lower bounds for \code{c(h2, intercept, eta, fold1, fold2, ..., foldn)}. Default is \code{NULL}, which means \code{c(0, 0, -2, 0, 0, ..., 0)}. When \code{fix.h2}, \code{fix.intercept} or \code{fix.eta} is not \code{NULL}, the lower bounds for \code{h2}, \code{intercept} or \code{eta} will be ignored.
#' @param upr Upper bounds for \code{c(h2, intercept, eta, fold1, fold2, ..., foldn)}, default is \code{NULL}, which means \code{c(1, 5, 1, M/Md_1, M/Md_2, ..., M/Md_n)}, where \code{Md_i} is the sum of the weights of the i-th annotation. When \code{fix.h2}, \code{fix.intercept} or \code{fix.eta} is not \code{NULL}, the upper bounds for \code{h2}, \code{intercept} or \code{eta} will be ignored.
#' @param pattern Chromosome and picece pattern of LD files, default is \code{".*chr(\\\\d{1,2})\\\\.(\\\\d{1,2})[_\\\\.].*"}.
#' @param par.h2 Whether to estimate the partitioned heritability, default \code{par.h2 = FALSE}.
#' @param nthreads Number of threads to use for matrix operations, default \code{nthreads = 1}. The default value is suitable for most cases, do not change it unless you are sure about the performance.
#' @return  A data.frame is returned with:
#' \itemize{
#' \item{item } The name of the parameter.
#' \item{estimation } The estimated value of the parameter.
#' \item{se } The standard error of the parameter.
#' \item{p } The p-value of the parameter.
#' \item{note } The note of the parameter.
#' }
#' @export
#'
#' @importFrom dplyr filter distinct
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads

HDL.S.joint <- function(
  Dr.path, gwas.df, LD.path,
  output.file = NULL, annots = NULL, mc.cores = 1, stepwise = FALSE,
  fill.missing.N = c("none", "min", "max", "median", "mean"),
  lim = exp(-18), lam.cut = NULL, verbose = FALSE,
  fix.h2 = NULL, fix.intercept = NULL, fix.eta = NULL, maxit=1000,
  pgtol=1e-3, start.v = NULL, lwr=NULL, upr=NULL,
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*", par.h2=FALSE, nthreads = 1){

  if(!is.null(output.file)){
    log.file <- paste0(output.file, ".log")
    if(file.exists(log.file)) unlink(log.file)
  }else{
    log.file <- ""
  }

  HDL.S:::log.msg("Starting HDL.S analysis...\n", log.file)
  gwas.df <- HDL.S:::format.gwas(gwas.df, LD.path, fill.missing.N, log.file, pattern)
  N <- median(gwas.df$N, na.rm=TRUE)
  z <- gwas.df$Z
  z <- matrix(z, ncol=1)
  rownames(z) <- gwas.df$SNP

  ref.data <- HDL.S:::load.Dr(Dr.path, annots, lam.cut, log.file, pattern, mc.cores)
  zr.data <- HDL.S:::HDL.S.reduct.dim(
    LD.path, z=z, D=NULL, lam.cut=lam.cut,
    nthreads=nthreads, pattern=pattern
  )
  for(i in seq_along(ref.data)){
    ref.data[[i]]$zr <- zr.data[[i]]$zr
  }

  HDL.S:::log.msg("Starting optimization...\n", log.file)
  res <- HDL.S:::HDL.S.optim(
    ref.data=ref.data, N=N, start.v=start.v,
    output.file=output.file, log.file=log.file, stepwise=stepwise,
    fix.h2=fix.h2, fix.intercept=fix.intercept, fix.eta=fix.eta,
    lim=lim, verbose=verbose, lwr=lwr, upr=upr,
    maxit=maxit, pgtol=pgtol, par.h2=par.h2,
    mc.cores=mc.cores, nthreads=nthreads
  )
  return(res)
}