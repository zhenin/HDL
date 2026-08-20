#' Perform heritability enrichment analysis of trait based on GWAS summary statistics.
#'
#' @param D A single-column matrix of annotation weights with rownames of SNP IDs and colname specifying the annotation name. Must have exactly one column.
#' @param gwas.df A data frame including GWAS summary statistics of genetic variants for a trait.
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS.
#' @param LD.path Path to the directory where linkage disequilibrium (LD) information is stored, compatible with HDL LD reference panels.
#' @param output.file Where the log and results should be written. If you do not specify a file, the log will be printed on the console.
#' @param mc.cores Number of cores to use for parallelization, default \code{mc.cores = 1}.
#' @param stepwise Whether to estimate enrichment fold by estimating heritability, intercept and eta first, default \code{stepwise = FALSE}. If \code{fix.h2} and \code{fix.intercept} are not \code{NULL}, \code{stepwise} will be overridden.
#' @param fill.missing.N If the sample size is missing in the GWAS summary statistics, \code{fill.missing.N} is used to fill the missing values.
#' @param lim Tolerance limitation to ensure the positive-definiteness of covariance matrices, default \code{lim = exp(-18)}.
#' @param lam.cut Eigenvalue cutoff for LD matrices, default \code{lam.cut = NULL}, which means no cutoff. For analyses with a limited number of traits and annotations, a lower cutoff (such as 0.1, or even not using a cutoff at all) is recommended. For large-scale analyses, a higher cutoff (such as 1)  is recommended, to yield fast computation.
#' @param Dr.path Path to the directory where the Dr matrices are stored, default \code{Dr.path = "./Dr"}.
#' @param overwrite Whether to overwrite the existing Dr matrices, default \code{overwrite = FALSE}.
#' @param verbose Whether to print the log on the console, default \code{verbose = FALSE}.
#' @param fix.h2 Whether to fix the heritability to \code{fix.h2} or estimate the heritability, default \code{fix.h2 = NULL}, which means estimate the heritability.
#' @param fix.intercept Whether to fix the intercept to \code{fix.intercept} or estimate the intercept, default \code{fix.intercept = NULL}, which means estimate the intercept.
#' @param fix.eta Whether to fix the eta to \code{fix.eta} or estimate the eta, default \code{fix.eta = NULL}, which means estimate the eta. If \code{fix.eta} is not \code{NULL}, the eta will be fixed to \code{fix.eta}.
#' @param maxit Maximum number of iterations, default \code{maxit = 1000}.
#' @param pgtol Tolerance for convergence, default \code{pgtol = 1e-3}.
#' @param start.v A vector of starting values \code{c(h2, intercept, eta, fold)} for optimization. Default is \code{NULL}, which means \code{c(0.1, 1, 0, 1)}. If \code{fix.h2}, \code{fix.intercept} or \code{fix.eta} is not \code{NULL}, the starting values for \code{h2}, \code{intercept} or \code{eta} will be ignored.
#' @param lwr Lower bounds for \code{c(h2, intercept, eta, fold)}. Default is \code{NULL}, which means \code{c(0, 0, -2, 0)}. When \code{fix.h2}, \code{fix.intercept} or \code{fix.eta} is not \code{NULL}, the lower bounds for \code{h2}, \code{intercept} or \code{eta} will be ignored.
#' @param upr Upper bounds for \code{c(h2, intercept, eta, fold)}, default is \code{NULL}, which means \code{c(1, 5, 1, M/Md)}, where \code{Md} is the sum of the weights of the annotation. When \code{fix.h2}, \code{fix.intercept} or \code{fix.eta} is not \code{NULL}, the upper bounds for \code{h2}, \code{intercept} or \code{eta} will be ignored.
#' @param mode If \code{mode = "disk"}, \code{Dr} is stored to disk (path returned only). If \code{mode = "memory"}, \code{Dr} is loaded to memory (matrix returned). Default is \code{mode = "disk"}.
#' @param pattern Chromosome and picece pattern of LD files, default is \code{".*chr(\\\\d{1,2})\\\\.(\\\\d{1,2})[_\\\\.].*"}.
#' @param norm.method The normalization method, either \code{"minmax"} (default), \code{"scaled"} or \code{"none"}. If \code{"minmax"}, the annotation weight vector \code{D} is normalized to [0, 1]. If \code{"scaled"}, the sum of normalized vector \code{D} is scaled to the number of annotated SNPs. If \code{"none"}, the annotation weight vector \code{D} is not normalized.
#' @param par.h2 Whether to estimate the partitioned heritability, default \code{par.h2 = FALSE}.
#' @param nthreads Number of threads to use for matrix operations, default \code{nthreads = 1}. The default value is suitable for most cases, do not change it unless you are sure about the performance.
#' @note Users can download the precomputed eigenvalues and eigenvectors of LD correlation matrices for European ancestry population. The download link can be found at https://github.com/zhenin/HDL/wiki/Reference-panels
#' These are the LD matrices and their eigen-decomposition from 335,265 genomic British UK Biobank individuals. Three sets of reference panel are provided:
#' 1) 1,029,876 QCed UK Biobank imputed HapMap3 SNPs. The size is about 33 GB after unzipping. Although it takes more time, using the imputed panel provides more accurate estimates of genetic correlations.
#' Therefore if the GWAS includes most of the HapMap3 SNPs, then we recommend using the imputed reference panel.
#' 2) 769,306 QCed UK Biobank imputed HapMap2 SNPs. The size is about 18 GB after unzipping.If one of your GWAS includes most of the HapMap 2 SNPs, but many SNPs (more than 1%) in the above HapMap 3 reference panel are absent,
#' then this HapMap2 panel is more proper to be used for HDL.
#' 3) 307,519 QCed UK Biobank Axiom Array SNPs. The size is about 7.5 GB after unzipping.
#'
#' @return  A data.frame is returned with:
#' \itemize{
#' \item{item } The name of the parameter.
#' \item{estimation } The estimated value of the parameter.
#' \item{se } The standard error of the parameter.
#' \item{p } The p-value of the parameter.
#' \item{note } The note of the parameter.
#' }
#' @author Ao Lan, Xia Shen
#'
#' @references
#' Lan A and Shen X (2024). Modeling the Genetic Architecture of Complex Traits via Stratified High-Definition Likelihood.
#'
#' @importFrom dplyr filter distinct
#' @importFrom parallel mclapply
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads
#'
#' @examples
#'\dontrun{
#' ## The GWAS summary statistics for birth weight loaded from HDL package.
#' data(gwas1.example, package='HDL')
#' M <- nrow(gwas1.example)
#' set.seed(1234)
#' D <- matrix(rbinom(M, 1, 0.01), ncol=1) # random D vector
#' row.names(D) <- gwas1.example$SNP
#' colnames(D) <- 'testAnno'
#'
#' ## The path to the directory where linkage disequilibrium (LD) information is stored.
#' LD.path <- "path/to/UKB_array_SVD_eigen90_extraction"
#' 
#' ## To speed up the test, we set a large lam.cut value.
#' res.HDL.S <- HDL.S(D, gwas1.example, LD.path, mc.cores=4, stepwise=TRUE, lam.cut=10, Dr.path=NULL, mode="memory")
#' print(as.data.frame(res.HDL.S))
#' }
#' @export
#'

HDL.S <- function(
  D, gwas.df, LD.path, output.file = NULL, mc.cores = 1, stepwise = FALSE,
  fill.missing.N = c("none", "min", "max", "median", "mean"),
  lim = exp(-18), lam.cut = NULL, Dr.path = "./Dr", overwrite = FALSE, verbose = FALSE,
  fix.h2 = NULL, fix.intercept = NULL, fix.eta = NULL, maxit=1000,
  pgtol=1e-3, start.v = NULL, lwr=NULL, upr=NULL, mode=c("disk", "memory"),
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*", norm.method=c("minmax", "scaled", "none"),
  par.h2=FALSE, nthreads = 1){

  mode <- match.arg(mode)
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

  HDL.S:::log.msg("Transfoming z (D) to zr (Dr)...\n", log.file)
  ref.data <- HDL.S:::HDL.S.reduct.dim(LD.path, z=z, D=D, lam.cut=lam.cut,
    Dr.path=Dr.path, overwrite=overwrite, mode=mode,
    mc.cores=mc.cores, pattern=pattern,
    norm.method=norm.method, log.file=log.file, nthreads=nthreads
  )
  HDL.S:::log.msg(paste0("Dr saved to ", Dr.path, "\n"), log.file)

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