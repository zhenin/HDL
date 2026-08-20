#' Reduce the dimension of covariance matrix by converting z (D) to zr (Dr).
#'
#' @param LD.path Path to the \code{.rda} file where the Eigen decomposition of LD matrix is stored.
#' @param D A matrix of annotation weights with rownames of SNP IDs and colnames specifying the annotation names.
#' @param z A matrix of Z-scores with rownames of SNP IDs. \bold{Supporting multiple columns for multiple traits.}
#' @param lam.cut Eigenvalue cutoff for LD matrices, default \code{lam.cut = NULL}, which means no cutoff. For analyses with a limited number of traits and annotations, a lower cutoff (such as 0.1, or even not using a cutoff at all) is recommended. For large-scale analyses, a higher cutoff (such as 1)  is recommended, to yield fast computation.
#' @param Dr.path Path to the directory where the Dr matrices are stored, default Dr.path = NULL, which means do not store Dr to disk.
#' @param overwrite Whether to overwrite the existing Dr matrices, default overwrite = FALSE.
#' @param mode Whether to store Dr to disk or memory, default \code{mode = "disk"}. If \code{mode = "disk"}, \code{Dr} is stored to disk (path returned only) and lam are not returned. If \code{mode = "memory"}, \code{Dr} and \code{lam} are returned.
#' @param mc.cores Number of cores to use for parallelization, default \code{mc.cores = 1}.
#' @param pattern Chromosome and picece pattern of LD files, default is \code{".*chr(\\\\d{1,2})\\\\.(\\\\d{1,2})[_\\\\.].*"}.
#' @param norm.method The normalization method, either \code{"minmax"} (default), \code{"scaled"} or \code{"none"}. If \code{"minmax"}, the annotation weight vector \code{D} is normalized to [0, 1]. If \code{"scaled"}, the sum of normalized vector \code{D} is scaled to the number of annotated SNPs. If \code{"none"}, the annotation weight vector \code{D} is not normalized.
#' @param log.file Where the log should be written. If you do not specify a file, the log will be printed on the console.
#' @param nthreads Number of threads to use for matrix operations, default \code{nthreads = 1}. The default value is suitable for most cases, do not change it unless you are sure about the performance.
#' @return A list is returned with:
#' \itemize{
#' \item{Dr }{The reduct RDR matrix.}
#' \item{zr }{The reduct z-score vector.}
#' \item{lam }{The eigenvalues of LD matrix.}
#' \item{Md }{The sum of annotation weights of SNPs in given LD matrix.}
#' \item{M }{The number of SNPs in given LD matrix.}
#' }
#' @export
#' 
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads

HDL.S.reduct.dim <- function(LD.path, z=NULL, D=NULL, lam.cut=NULL,
  Dr.path=NULL, overwrite=FALSE, mode=c("disk", "memory"),
  mc.cores=1, pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*",
  norm.method=c("minmax", "scaled", "none"), log.file="", nthreads=1){
  # use RhpcBLASctl to control the number of thread
  RhpcBLASctl::blas_set_num_threads(nthreads)
  RhpcBLASctl::omp_set_num_threads(nthreads)

  if(!is.null(D) && is.null(Dr.path)){
    Dr.path <- "./HDL.S_Dr"
    msg <- sprintf(
      "D is not NULL and Dr.path is not provided. The default path %s will be used.\n",
      Dr.path
    )
    HDL.S:::log.msg(msg, log.file)
  }

  if(!is.null(Dr.path) && Dr.path!="" && !dir.exists(Dr.path))
    dir.create(Dr.path, recursive = TRUE)
  mode <- match.arg(mode)
  if(is.null(Dr.path)) mode <- "memory"

  LD.files <- HDL.S:::list.LD.ref.files(LD.path, suffix=".rda",
    pattern=pattern, log.file=log.file)

  if(!is.null(D)) D <- HDL.S:::normD(D, LD.path, log.file, norm.method, pattern)

  ref.data <- parallel::mclapply(
    LD.files, HDL.S:::trans.zD, Dr.path=Dr.path, z=z, D=D,
    lam.cut=lam.cut, mode=mode,
    overwrite=overwrite, log.file=log.file, nthreads=nthreads,
    mc.cores=mc.cores
  )

  return(ref.data)
}