#' @title log messages and errors
#' @noRd
#' @keywords internal
log.msg <- function(message, log.file = "", type="message"){
  if(is.null(message)) return(NULL)
  if(log.file != ""){
    cat(message, file = log.file, append = T)
  }
  if(type == "warning"){
    warning(message)
  }else if(type == "error"){
    stop(message)
  }else if(log.file == "" && type == "message"){
    cat(message)
  }
}

#' @title Format GWAS data
#' @noRd
#' @keywords internal
format.gwas <- function(gwas.df, LD.path,
  fill.missing.N = c("none", "min", "max", "median", "mean"), log.file = "",
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"){
  fill.missing.N <- match.arg(fill.missing.N)
  bim.files <- HDL.S:::list.LD.ref.files(LD.path, suffix='\\.bim',
    pattern=pattern, log.file=log.file)
  bim <- do.call(rbind, lapply(
    bim.files, read.table, header=F,
    col.names=c('CHR', 'SNP', 'CM', 'POS', 'A1', 'A2')
  ))
  Mref <- nrow(bim)

  gwas.df <- filter(gwas.df, SNP %in% bim$SNP)

  gwas.df$A1 <- toupper(as.character(gwas.df$A1))
  gwas.df$A2 <- toupper(as.character(gwas.df$A2))

  Z.found <- "Z" %in% colnames(gwas.df)
  b.found <- "b" %in% colnames(gwas.df)
  se.found <- "se" %in% colnames(gwas.df)
  if(!Z.found && !b.found && !se.found){
    err.msg <- "Z, b or se is missing in GWAS. Please check.\n"
    HDL.S:::log.msg(err.msg, log.file, type="error")
  }
  if(!Z.found && b.found && se.found){
    if(abs(median(gwas.df$b) - 1) < 0.1){
      warn.msg <- "Taking log(b) in GWAS because b is likely to be OR in stead of log(OR). \n"
      HDL.S:::log.msg(warn.msg, log.file, type="warning")
      gwas.df$Z <- log(gwas.df$b) / gwas.df$se
    } else{
      gwas.df$Z <- gwas.df$b / gwas.df$se
    }
  }

  gwas.df <- filter(gwas.df, !is.na(Z), is.finite(Z))
  gwas.df$N[is.na(gwas.df$N)] <- switch(
    fill.missing.N,
    "min" = min(gwas.df$N, na.rm = T),
    "max" = max(gwas.df$N, na.rm = T),
    "median" = median(gwas.df$N, na.rm = T),
    "mean" = mean(gwas.df$N, na.rm = T),
    NA
  )
  gwas.df <- filter(gwas.df, !is.na(N))
  gwas.M <- length(unique(gwas.df$SNP))
  if(gwas.M < Mref*0.99){
    warn.msg <- "More than 1%% SNPs in reference panel are missed in the GWAS. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
    HDL.S:::log.msg(warn.msg, log.file, type="warning")
  }

  HDL.S:::log.msg(
    sprintf(
      "%d out of %d (%.2f%%) SNPs in reference panel are available in the GWAS.\n",
      gwas.M, Mref, 100*gwas.M/Mref
    ),
    log.file
  )

  ## remove duplicates
  gwas.df <- distinct(gwas.df, SNP, .keep_all = TRUE)

  ## match direction of z-scores
  ref.A1 <- bim$A1
  names(ref.A1) <- bim$SNP
  ref.A1 <- ref.A1[gwas.df$SNP]
  idx.direc <- ifelse(ref.A1 == gwas.df$A1, 1, -1)
  gwas.df$Z <- gwas.df$Z * idx.direc
  return(gwas.df)
}

#' @title List LD reference files
#' @noRd
#' @keywords internal
list.LD.ref.files <- function(LD.path, suffix = ".rda", full.names = TRUE,
                              pattern = ".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*",
                              log.file=""){
  all.segfiles <- list.files(LD.path, pattern = suffix, full.names = F)
  segfiles <- all.segfiles[grepl(pattern = pattern, all.segfiles)]
  chrom <- gsub(pattern, "\\1", segfiles)
  piece <- gsub(pattern, "\\2", segfiles)
  if(full.names) segfiles <- paste0(LD.path, "/", segfiles)
  if(length(chrom)==0){
    
    if(length(all.segfiles) > 0){
      msg <- 'No valid LD reference files found matched the pattern. All files with given suffix will be used. \n'
      HDL.S:::log.msg(msg, log.file, type="warning")
      if(full.names) return(paste0(LD.path, "/", all.segfiles))
    }
    msg <- 'No valid LD reference files found. Please check the suffix and pattern. \n'
    HDL.S:::log.msg(msg, log.file, type="error")
    stop(-1)
  }
  chrom <- as.numeric(chrom)
  piece <- as.numeric(piece)
  idx <- order(chrom, piece)
  segfiles <- segfiles[idx]
  return(segfiles)
}

#' @title log likelihood function without enrichment fold parameter
#' @noRd
#' @keywords internal
log.lik.HDL <- function(param, lam, zr, N, M, lim=exp(-18),
  fix.h2=NULL, fix.intercept=NULL, fix.eta=NULL){
  pp <- rep(NA, 3)
  if(!is.null(fix.h2)) pp[1] <- fix.h2
  if(!is.null(fix.intercept)) pp[2] <- fix.intercept
  if(!is.null(fix.eta)) pp[3] <- fix.eta
  pp[is.na(pp)] <- param
  h2 <- pp[1]
  intercept <- pp[2]
  eta <- pp[3]
  alpha <- eta**2 * N * h2 / M - eta +
    (intercept + eta - 2 * eta * (intercept + eta) * N * h2 / M) * lam +
    (intercept + eta)**2 * N * h2 / M * lam**2
  alpha <- pmax(alpha, lim)
  logdet <- sum(log(alpha))
  quadra <- sum(zr**2/alpha)
  lnL <- -0.5 * (logdet + quadra)
  return(lnL)
}

#' @title log likelihood function for whole genome
#' @noRd
#' @keywords internal
log.lik.wg <- function(param, ref.data, Md, M, N,
  log.file="", fix.h2=NULL, fix.intercept=NULL, fix.eta=NULL,
  verbose=FALSE, lim=exp(-18), mc.cores=1, nthreads=1){
  t0 <- Sys.time()

  folds <- param[seq(length(param) - length(Md) + 1, length(param))]
  params <- param[seq_len(length(param) - length(Md))]
  pp <- rep(NA, 3)
  if(!is.null(fix.h2)) pp[1] <- fix.h2
  if(!is.null(fix.intercept)) pp[2] <- fix.intercept
  if(!is.null(fix.eta)) pp[3] <- fix.eta
  pp[is.na(pp)] <- params
  h2 <- pp[1]
  intercept <- pp[2]
  eta <- pp[3]

  per.h2d <- (folds - 1) / (M - Md) * h2
  per.h20 <- (h2 - sum(Md * per.h2d)) / M

  if(mc.cores>1){
    lnL <- parallel::mclapply(
      ref.data, HDL.S:::log.lik,
      per.h20=per.h20, per.h2d=per.h2d, intercept=intercept, eta=eta,
      N=N, lim=lim, nthreads=nthreads, mc.cores=mc.cores
    )
  }else{
    lnL <- lapply(
      ref.data, HDL.S:::log.lik,
      per.h20=per.h20, per.h2d=per.h2d, intercept=intercept, eta=eta,
      N=N, lim=lim, nthreads=nthreads
    )
  }

  lnL <- sum(as.numeric(unlist(lnL)))
  time <- as.numeric(Sys.time() - t0, units="secs")
  if(verbose){
    folds <- paste0(round(folds, 3), collapse=", ")
    msg <- sprintf(
      "fold(s): %s h2: %.3f intercept: %.3f eta: %.3f lnL: %.3f time: %.3f \n",
      folds, h2, intercept, eta, lnL, time
    )
    HDL.S:::log.msg(msg, log.file)
  }
  return(lnL)
}

#' @title D normalization
#' @noRd
#' @keywords internal
normD <- function(
  D, LD.path, log.file="", norm.method=c("minmax", "scaled", "none"),
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"){

  method <- match.arg(norm.method)
  msg <- paste0('Annotation name(s): [', paste0(colnames(D), collapse = ', '), '].\n')
  HDL.S:::log.msg(msg, log.file)
  bim.files <- HDL.S:::list.LD.ref.files(LD.path, suffix='\\.bim',
    pattern=pattern, log.file=log.file)
  bim <- do.call(rbind, lapply(
    bim.files, read.table, header=F,
    col.names=c('CHR', 'SNP', 'CM', 'POS', 'A1', 'A2')
  ))
  ref.snps <- bim$SNP
  M <- nrow(bim)
  dup.snps <- row.names(D)[duplicated(row.names(D))]
  D <- D[setdiff(row.names(D), dup.snps), , drop=F] ## remove duplicates

  int.snps <- intersect(ref.snps, row.names(D))
  
  
  D <- D[int.snps, , drop=F]
  dD <- matrix(0, nrow=M, ncol=ncol(D), dimnames=(list(ref.snps, colnames(D))))
  dD[int.snps, ] <- D[int.snps, ]
  dD[is.na(dD)] <- 0
  D <- dD

  msg.prefix <- sprintf('Applied `%s` weight nomalization', method)
  minv <- apply(D, 2, min, na.rm=T)
  if(method == "minmax"){
    maxv <- apply(D, 2, max, na.rm=T)
    D <- t((t(D) - minv) / (maxv - minv))
    Md <- colSums(D, na.rm=T)
  }else if(method == "scaled"){
    minv[minv > 0] <- 0
    D <- t(t(D) - minv)
    Md <- colSums(D > 0, na.rm=T)
    idx <- Md == M
    if(sum(idx) >= 1){
      minv <- apply(D, 2, min, na.rm=T)
      D <- t(t(D) - minv)
      Md <- colSums(D > 0, na.rm=T)
    }
    D <- t(Md / colSums(D, na.rm=T) * t(D))
  }else if(method=="none"){
    if (sum(minv < 0) >= 1){
      neg.anns <- paste0(names(minv)[minv < 0], collapse=',')
      warn.msg <- sprintf("The annotation (%s) weights contain negative values, which may cause error (or bias) in the estimation.\n", neg.anns)
      HDL.S:::log.msg(warn.msg, log.file, type="warning")
    }
    msg.prefix <- 'No normalization applied'
  }else{
    err.msg <- "Unknown normalization method. Please choose one of 'minmax', 'scaled', 'none'.\n"
    HDL.S:::log.msg(err.msg, log.file, type="error")
  }

  Ma <- colSums(D != 0, na.rm=T)
  Md <- colSums(D, na.rm=T)
  upr <- paste0(sapply(M / Md, function(x) sprintf('%.3f', x)), collapse=', ')
  cnt <- paste0(sapply(seq_along(Ma),
    function(x) sprintf('%d (%.3f%%)', Ma[x], Ma[x]/M*100)), collapse=', ')
  msg <- sprintf(
    "%s on [%s] annotated variants. The theoretical upper boundary for enrichment fold is M / Md = [%s].\n",
    msg.prefix, cnt, upr
  )
  HDL.S:::log.msg(msg, log.file)

  if(sum(Md/M > 0.9) >= 1){
    warn.annots <- paste0(names(Md)[Md / M > 0.9], collapse=', ')
    warn.msg <- sprintf("Upper boundary for the annotation [%s] is less than or very close to 1, please consider reducing the annotation density or using other normalization methods.\n", warn.annots)
    HDL.S:::log.msg(warn.msg, log.file, type="warning")
  }

  return(dD)
}

#' @title log likelihood function for joint optimization
#' @noRd
#' @keywords internal
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom RhpcBLASctl omp_set_num_threads
log.lik <- function(refd, per.h20, per.h2d, intercept, eta,
  N, lim = exp(-18), nthreads = 1){
  RhpcBLASctl::blas_set_num_threads(nthreads)
  RhpcBLASctl::omp_set_num_threads(nthreads)
  zr <- refd$zr
  lam <- refd$lam
  
  Dr <- HDL.S:::sum.Dr(refd$Dr, per.h2d, length(lam))

  alpha <- eta**2 * N * per.h20 - eta +
    (intercept + eta - 2 * eta * (intercept + eta) * N * per.h20) * lam +
    (intercept + eta)**2 * N * per.h20 * lam**2
  alpha <- pmax(alpha, lim)

  if(length(Dr) == 1){
    Dr <- as.vector(Dr)
    sigma <- alpha + N * Dr
    logdet <- log(sigma)
    quadra <- zr**2/sigma
    lnL <- -0.5 * (logdet + quadra)
    return(lnL)
  }

  if(eta == 0){
    sigma <- (intercept**2 * N) * Dr
  }else{
    ss <- -eta * (intercept + eta) * N / lam
    # sigma <- diag(ss + (intercept + eta)**2 * N) %*% Dr +
    #   Dr %*% diag(ss) + eta**2 * N * (diag(1/lam) %*% Dr %*% diag(1/lam))
    cc <- (intercept + eta)**2 * N
    inv.lam <- 1 / lam
    sigma <- Dr * (ss + cc) + t(t(Dr) * ss) + 
      eta**2 * N * Dr * outer(inv.lam, inv.lam)
  }
  diag(sigma) <- diag(sigma) + alpha
  
  L <- tryCatch({
    chol(sigma)
  }, error = function(e){
    NA
  })
  if(length(L)==1 && is.na(L)) return(-1e18)
  logdet <- tryCatch({
    2 * sum(log(diag(L)))
  }, error = function(e){
    idx <- diag(L) > 0
    if(sum(idx)==0) return(NA)
    return(2 * sum(log(diag(L)[idx])))
  })
  quadra <- crossprod(zr,
    forwardsolve(L, backsolve(L, zr, transpose=T), upper.tri=T))
  logdet <- as.vector(logdet)
  quadra <- as.vector(quadra)
  lnL <- -0.5 * (logdet + quadra)
  if(is.na(lnL)) return(-1e18)
  return(lnL)
}

#' @title partitioned heritability estimation
#' @noRd
#' @keywords internal
est.par.h2 <- function(h2, h2.se, folds, folds.se, cov.h2.folds, Mds, M){
  if(length(folds) != length(Mds) || length(folds) != length(cov.h2.folds)){
    HDL.S:::log.msg('The length of folds, Mds and cov.h2.folds should be the same.\n', type='error')
  }
  a <- Mds / M

  ph2 <- a * folds * h2
  d.ph2.h2 <- a * folds
  d.ph2.folds <- a * h2

  v.h2 <- h2.se**2
  v.folds <- folds.se**2
  v.ph2 <- d.ph2.h2**2 * v.h2 + 
    d.ph2.folds**2 * v.folds +
    2 * d.ph2.h2 * d.ph2.folds * cov.h2.folds
  v.ph2[v.ph2 < 0] <- NA
  ph2.se <- sqrt(v.ph2)

  return(list(ph2=ph2, ph2.se=ph2.se))
}

#' @title sum Dr
#' @noRd
#' @keywords internal
sum.Dr <- function(Dr, per.h2d, r){
  idx <- seq_len(r)
  if(length(per.h2d)==1){
    cur.Dr <- Dr[[1]]
    if(length(cur.Dr)==1 && is.character(cur.Dr) && file.exists(cur.Dr)){
      e1 <- new.env()
      load(cur.Dr, envir = e1)
      cur.Dr <- e1$Dr
    }
    return(cur.Dr[idx, idx, drop=F] * per.h2d)
  }
  sDr <- NULL
  for(i in seq_along(per.h2d)){
    cur.Dr <- Dr[[i]]
    if(length(cur.Dr)==1 && is.character(cur.Dr) && file.exists(cur.Dr)){
      e1 <- new.env()
      load(cur.Dr, envir = e1)
      cur.Dr <- e1$Dr
    }
    cur.Dr <- cur.Dr[idx, idx, drop=F] * per.h2d[i]
    if(is.null(sDr)){
      sDr <- cur.Dr
    }else{
      sDr <- sDr + cur.Dr
    }
  }
  return(sDr)
}

#' @title load ref.data from Dr.path
#' @noRd
#' @keywords internal
load.Dr <- function(Dr.path, kept.annots=NULL, lam.cut=NULL, log.file="",
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*", mc.cores=1, mode=c("disk", "memory")){
  if(!file.exists(Dr.path)){
    err.msg <- sprintf("The file %s does not exist.\n", Dr.path)
    HDL.S:::log.msg(err.msg, log.file, type="error")
  }
  mode <- match.arg(mode)
  Dr.path <- normalizePath(Dr.path)
  annots <- list.dirs(Dr.path, full.names=FALSE, recursive=FALSE)
  if(length(annots)==0 && is.null(kept.annots)){
    annots <- basename(Dr.path)
    Dr.path <- dirname(Dr.path)
  }else if(!is.null(kept.annots)){
    if(length(intersect(annots, kept.annots)) != length(kept.annots)){
      err.msg <- sprintf("The required annotations %s are not found in the Dr.path %s.\n",
        paste(setdiff(kept.annots, annots), collapse=", "), Dr.path)
      HDL.S:::log.msg(err.msg, log.file, type="error")
    }
    annots <- intersect(annots, kept.annots)
  }

  Dr.files <- HDL.S:::list.LD.ref.files(file.path(Dr.path, annots[1]), suffix='.rda',
    pattern=pattern, log.file=log.file)
  Dr.files <- basename(Dr.files)
  ref.data <- parallel::mclapply(
    Dr.files, function(Dr.file, annots, Dr.path, lam.cut){
      rda.files <- paste0(Dr.path, '/', annots, '/', Dr.file)
      e1 <- new.env()
      load(rda.files[1], envir = e1)
      lam <- e1$lam
      Mds <- e1$Md
      M <- e1$M

      if(!is.null(lam.cut)){
        idx <- lam > lam.cut
        if(sum(idx)==0) idx[1] <- TRUE
        lam <- lam[idx]
      }
      idx <- seq_len(length(lam))

      Dr <- as.list(paste0(Dr.path, '/', annots, '/', Dr.file))
      names(Dr) <- annots
      if(mode=='memory') Dr[[1]] <- e1$Dr[idx, idx, drop=F]
      if(length(rda.files) == 1){
        names(Mds) <- annots
        return(list(Dr=Dr, lam=lam, Md=Mds, M=M))
      }
      for(j in seq(2, length(rda.files))){
        e1 <- new.env()
        load(rda.files[j], envir = e1)
        Mds <- c(Mds, e1$Md)
        if(mode=='memory') Dr[[j]] <- e1$Dr[idx, idx, drop=F]
      }
      names(Mds) <- annots
      return(list(Dr=Dr, lam=lam, Md=Mds, M=M))
    }, annots=annots, Dr.path=Dr.path, lam.cut=lam.cut, mc.cores=mc.cores
  )
  return(ref.data)
}


#' @title subset D or z matrix
#' @noRd
#' @keywords internal
subset.zD <- function(V, zD){
  snps.ref <- row.names(V)
  sub.zD <- matrix(0, nrow=nrow(V), ncol=ncol(zD))
  row.names(sub.zD) <- snps.ref
  colnames(sub.zD) <- colnames(zD)
  int.snps <- intersect(snps.ref, row.names(zD))
  sub.zD[int.snps,] <- zD[int.snps,]
  sub.zD[is.na(sub.zD)] <- 0
  sub.zD[!is.finite(sub.zD)] <- 0
  return(sub.zD)
}


#' @title transfrom z and D to zr and Dr
#' @noRd
#' @keywords internal
trans.zD <- function(
  LD.file,
  Dr.path=NULL, z=NULL, D=NULL, lam.cut=NULL,
  mode=c("disk", "memory"), overwrite=FALSE, log.file="", nthreads=1){

  RhpcBLASctl::blas_set_num_threads(nthreads)
  RhpcBLASctl::omp_set_num_threads(nthreads)
  mode <- match.arg(mode)
  if(!is.null(D) && is.null(Dr.path)){
    Dr.path <- "./HDL.S_Dr"
    msg <- sprintf(
      "D is not NULL and Dr.path is not provided. The default path %s will be used.\n",
      Dr.path
    )
    HDL.S:::log.msg(msg, log.file)
  }

  load(LD.file) # V, lam, LDsc
  if(exists('U') && !exists('V')){
    V <- U
    rm('U')
  }
  if(!is.null(lam.cut)){
    idx <- lam > lam.cut
    if(sum(idx)==0) idx[1] <- TRUE
    V <- V[,idx]
    lam <- lam[idx]
  }
  snps.ref <- row.names(V)
  M <- length(snps.ref)

  zr <- NULL
  ## convert z to zr
  if(!is.null(z)) zr <- crossprod(V, HDL.S:::subset.zD(V, z))

  if(!is.null(D)){
    if(!dir.exists(Dr.path)) dir.create(Dr.path, recursive=TRUE)
    D <- HDL.S:::subset.zD(V, D) ## subset D to reference panel
    Mds <- colSums(D)
    names(Mds) <- colnames(D)
    M <- nrow(D)
    annots <- colnames(D)
    Dr.files <- paste0(Dr.path, '/', annots, '/', basename(LD.file))
    names(Dr.files) <- annots
    for(i in seq_along(annots)){
      annot <- annots[i]
      Dr.file <- Dr.files[i]
      if(!dir.exists(dirname(Dr.file))) dir.create(dirname(Dr.file), recursive=TRUE)
      if(!file.exists(Dr.file) || overwrite){
        ## save Dr to disk
        Dr <- diag(lam) %*% t(V) %*% diag(D[, i, drop=TRUE]) %*% V %*% diag(lam)
        Md <- Mds[annot]
        save(Dr, lam, Md, M, file=Dr.file)
      }
    }
  }

  if(is.null(D) && !is.null(Dr.path)){
    annots <- list.dirs(Dr.path, full.names=FALSE, recursive=FALSE)
    if(length(annots)==0){
      annots <- basename(Dr.path)
      Dr.path <- dirname(normalizePath(Dr.path))
    }
    Dr.files <- paste0(Dr.path, '/', annots, '/', basename(LD.file))
    names(Dr.files) <- annots
  }

  ref.d <- list(
    Dr=NULL,
    lam=lam,
    Md=NULL,
    M=M,
    zr=NULL
  )
  if(!is.null(zr)){
    ref.d$zr <- zr
  }
  if(!is.null(D) || !is.null(Dr.path)){
    ref.d$Dr <- as.list(Dr.files)

    Mds <- c()
    for(i in seq_len(length(ref.d$Dr))){
      e1 <- new.env()
      load(ref.d$Dr[[i]], envir = e1)
      if(mode=="memory") ref.d$Dr[[i]] <- e1$Dr
      Mds <- c(Mds, e1$Md)
    }
    names(Mds) <- names(ref.d$Dr)
    ref.d$Md <- Mds
    ref.d$M <- e1$M
  }

  return(ref.d)
}