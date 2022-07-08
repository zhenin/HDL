#' High-definition likelihood inference of genetic correlations (HDL)
#' 
#' The function returns the estimate and standard error of the genetic correlation between two traits based on GWAS summary statistics. 
#' 
#' @param gwas1.df A data frame including GWAS summary statistics of genetic variants for trait 1. 
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS. 
#' @param gwas2.df A data frame including GWAS summary statistics of genetic variants for trait 2. 
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS.
#' @param LD.path Path to the directory where linkage disequilibrium (LD) information is stored.
#' @param Nref Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 335265
#' @param N0 Number of individuals included in both cohorts. The estimated genetic correlation is usually robust against misspecified N0. 
#' If not given, the default value is set to the minimum sample size across all SNPs in cohort 1 and cohort 2.
#' @param output.file Where the log and results should be written. If you do not specify a file, the log will be printed on the console.
#' @param eigen.cut Which eigenvalues and eigenvectors in each LD score matrix should be used for HDL. 
#' Users are allowed to specify a numeric value between 0 and 1 for eigen.cut. For example, eigen.cut = 0.99 means using the leading eigenvalues explaining 99% of the variance
#' and their correspondent eigenvectors. If the default 'automatic' is used, the eigen.cut gives the most stable heritability estimates will be used. 
#' @param lim Tolerance limitation, default lim = exp(-18). 
#' @note Users can download the precomputed eigenvalues and eigenvectors of LD correlation matrices for European ancestry population. The download link can be found at https://github.com/zhenin/HDL/wiki/Reference-panels
#' These are the LD matrices and their eigen-decomposition from 335,265 genomic British UK Biobank individuals. Two sets of reference panel are provided: 
#' 1) 307,519 QCed UK Biobank Axiom Array SNPs. The size is about 7.5 GB after unzipping.
#' 2) 1,029,876 QCed UK Biobank imputed SNPs. The size is about 31 GB after unzipping. Although it takes more time, using the imputed panel provides more accurate estimates of genetic correlations. 
#' Therefore if the GWAS includes most of the HapMap3 SNPs, then we recommend using the imputed reference panel.
#' 
#' 
#' @return A list is returned with:
#' \itemize{
#' \item{rg }{The estimated genetic correlation.}
#' \item{rg.se }{The standard error of the estimated genetic correlation.}
#' \item{P }{P-value based on Wald test.}
#' \item{estimates.df }{A detailed matrix includes the estimates and standard errors of heritabilities, genetic covariance and genetic correlation.}
#' \item{eigen.use }{The eigen.cut used in computation.}
#' }
#' 
#' @author Zheng Ning, Xia Shen
#' 
#' @references 
#' Ning Z, Pawitan Y, Shen X (2020). High-definition likelihood inference of genetic correlations 
#' across human complex traits. \emph{Accepted by Nature Genetics, waiting for publishing}.
#' 
#' @seealso 
#' HDL tutorial: https://github.com/zhenin/HDL
#' @import dplyr
#' 
#' @examples 
#'\dontrun{
#' ## The GWAS summary statistics for birth weight 
#' data(gwas1.example)
#' 
#' ## The GWAS summary statistics for type 2 diabetes
#' data(gwas2.example)
#' 
#' ## The path to the directory where linkage disequilibrium (LD) information is stored.
#' LD.path <- "/Users/zhengning/Work/HDL/package/UKB_array_SVD_eigen90_extraction"
#' 
#' res.HDL <- HDL.rg(gwas1.example, gwas2.example, LD.path)
#' res.HDL
#' }
#' @export
#' 

HDL.rg <-
  function(gwas1.df, gwas2.df, LD.path, Nref = 335265, N0 = min(gwas1.df$N, gwas2.df$N), output.file = "", eigen.cut = "automatic", 
           jackknife.df = FALSE, intercept.output = FALSE, fill.missing.N = NULL, lim = exp(-18)){
    
    if(output.file != ""){
      if(file.exists(output.file) == T){
        system(paste0("rm ",output.file))
      }
      pkgDescription <- packageDescription("HDL")
      pkgVersion <- pkgDescription$Version
      pkgDate <- pkgDescription$Date
      pkgName <- pkgDescription$Package
      pkgTitle <- pkgDescription$Title
      pkgAuthor <- pkgDescription$Author
      pkgMaintainer <- pkgDescription$Maintainer
      cat(paste("\n", pkgName, ": ", pkgTitle, "\n", sep = ""), file = output.file, append = T)
      cat(paste("Version ", pkgVersion, " (", pkgDate, ") installed", "\n", sep = ""), file = output.file, append = T)
      cat(paste("Author: ", pkgAuthor, "\n", sep = ""), file = output.file, append = T)
      cat(paste("Maintainer: ", pkgMaintainer, "\n", sep = ""), file = output.file, append = T)
      cat('Tutorial: https://github.com/zhenin/HDL\n', file = output.file, append = T)
      cat('Use citation("HDL") to know how to cite this work.\n\n', file = output.file, append = T)
    }
    time.start <- date()
    cat("Analysis starts on",time.start,"\n")
    if(output.file != ""){
      cat("Analysis starts on",time.start,"\n", file = output.file, append = T)
    }
    
    if(eigen.cut != "automatic" && is.na(as.numeric(eigen.cut))){
      error.message <- "The input of eigen.cut has to be 'automatic' or a number between 0 and 1. \n"
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
    
    LD.files <- list.files(LD.path)
    
    if(any(grepl(x = LD.files, pattern = ".*_snp_counter.*"))){
      snp_counter_file <- LD.files[grep(x = LD.files, pattern = ".*_snp_counter.*")]
      snp_list_file <- LD.files[grep(x = LD.files, pattern = ".*_snp_list.*")]
      load(file=paste(LD.path, snp_counter_file, sep = "/"))
      load(file=paste(LD.path, snp_list_file, sep = "/"))
      if("nsnps.list.imputed" %in% ls()){
        snps.name.list <- snps.list.imputed.vector
        nsnps.list <- nsnps.list.imputed
      }
      if(is.null(names(nsnps.list))) names(nsnps.list) <- as.character(1:length(nsnps.list))
    } else{
      error.message <- "It seems this directory does not contain all files needed for HDL. Please check your LD.path again. The current version of HDL only support pre-computed LD reference panels."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
    
    gwas1.df <- gwas1.df %>% filter(SNP %in% snps.name.list)
    gwas2.df <- gwas2.df %>% filter(SNP %in% snps.name.list)
    
    gwas1.df$A1 <- toupper(as.character(gwas1.df$A1))
    gwas1.df$A2 <- toupper(as.character(gwas1.df$A2))
    gwas2.df$A1 <- toupper(as.character(gwas2.df$A1))
    gwas2.df$A2 <- toupper(as.character(gwas2.df$A2))
    
    
    if(!("Z" %in% colnames(gwas1.df))){
      if(("b" %in% colnames(gwas1.df)) && ("se" %in% colnames(gwas1.df))){
        if(abs(median(gwas1.df$b) - 1) < 0.1){
          cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n")
          if(output.file != ""){
            cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n", file = output.file, append = T)
          }
          gwas1.df$Z <- log(gwas1.df$b) / gwas1.df$se
        } else{
          gwas1.df$Z <- gwas1.df$b / gwas1.df$se
        }
      } else{
        error.message <- "Z is not available in GWAS 1, meanwhile either b or se is missing. Please check."
        if(output.file != ""){
          cat(error.message, file = output.file, append = T)
        }
        stop(error.message)
        
      }
    }
    
    if(!("Z" %in% colnames(gwas2.df))){
      if(("b" %in% colnames(gwas2.df)) && ("se" %in% colnames(gwas2.df))){
        if(abs(median(gwas2.df$b) - 1) < 0.1){
          cat("Taking log(b) in GWAS 2 because b is likely to be OR in stead of log(OR). \n")
          if(output.file != ""){
            cat("Taking log(b) in GWAS 2 because b is likely to be OR in stead of log(OR). \n", file = output.file, append = T)
          }
          gwas2.df$Z <- log(gwas2.df$b) / gwas2.df$se
        } else{
          gwas2.df$Z <- gwas2.df$b / gwas2.df$se
        }
      } else{
        error.message <- "Z is not available in GWAS 2, meanwhile either b or se is missing. Please check."
        if(output.file != ""){
          cat(error.message, file = output.file, append = T)
        }
        stop(error.message)
        
      }
    }
    
    k1.0 <- length(unique(gwas1.df$SNP))
    k2.0 <- length(unique(gwas2.df$SNP))
    
    if(is.null(fill.missing.N)){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z), !is.na(N))
      gwas2.df <- gwas2.df %>% filter(!is.na(Z), !is.na(N))
    } else if(fill.missing.N == "min"){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z))
      gwas1.df$N[is.na(gwas1.df$N)] <- min(gwas1.df$N, na.rm = T)
      
      gwas2.df <- gwas2.df %>% filter(!is.na(Z))
      gwas2.df$N[is.na(gwas2.df$N)] <- min(gwas2.df$N, na.rm = T)
    } else if(fill.missing.N == "max"){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z))
      gwas1.df$N[is.na(gwas1.df$N)] <- max(gwas1.df$N, na.rm = T)
      
      gwas2.df <- gwas2.df %>% filter(!is.na(Z))
      gwas2.df$N[is.na(gwas2.df$N)] <- max(gwas2.df$N, na.rm = T)
    } else if(fill.missing.N == "median"){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z))
      gwas1.df$N[is.na(gwas1.df$N)] <- median(gwas1.df$N, na.rm = T)
      
      gwas2.df <- gwas2.df %>% filter(!is.na(Z))
      gwas2.df$N[is.na(gwas2.df$N)] <- median(gwas2.df$N, na.rm = T)
    } else{
      error.message <- "If given, the argument fill.missing.N can only be one of below: 'min', 'max', 'median'."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
    
    k1 <- length(unique(gwas1.df$SNP))
    k2 <- length(unique(gwas2.df$SNP))
    k1.percent <- paste("(",round(100*k1 / length(snps.name.list), 2), "%)", sep="") 
    k2.percent <- paste("(",round(100*k2 / length(snps.name.list), 2), "%)", sep="") 
    
    cat(k1.0 - k1, "SNPs were removed in GWAS 1 due to missing N or missing test statistic."," \n")
    cat(k2.0 - k2, "SNPs were removed in GWAS 2 due to missing N or missing test statistic."," \n")
    if(output.file != ""){
      cat(k1.0 - k1, " SNPs were removed in GWAS 1 due to missing N or missing test statistic."," \n", file = output.file, append = T)
      cat(k2.0 - k2, " SNPs were removed in GWAS 2 due to missing N or missing test statistic."," \n", file = output.file, append = T)
    }
    
    cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in GWAS 1."," \n")
    cat(k2, "out of", length(snps.name.list), k2.percent, "SNPs in reference panel are available in GWAS 2."," \n")
    if(output.file != ""){
      cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in GWAS 1."," \n", file = output.file, append = T)
      cat(k2, "out of", length(snps.name.list), k2.percent, "SNPs in reference panel are available in GWAS 2."," \n", file = output.file, append = T)
    }
    if(k1 < length(snps.name.list)*0.99){
      error.message <- "Warning: More than 1% SNPs in reference panel are missed in GWAS 1. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      cat(error.message)
    }
    if(k2 < length(snps.name.list)*0.99){
      error.message <- "Warning: More than 1% SNPs in reference panel are missed in GWAS 2. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      cat(error.message)
    }
    
    
    
    ## stats
    N1 <- median(gwas1.df$N)
    N2 <- median(gwas2.df$N)
    N <- sqrt(N1)*sqrt(N2)
    
    ## samples for phenotypes
    p1 <- N0/N1
    p2 <- N0/N2
    
    rho12 <- suppressWarnings(inner_join(gwas1.df %>% select(SNP, Z), gwas2.df %>% select(SNP, Z), by = "SNP") %>%
                                summarise(x=cor(Z.x, Z.y, use = "complete.obs")) %>% unlist)
    
    bstar1.v <- bstar2.v <- lam.v <- list()
    HDL11.df <- HDL12.df <- HDL22.df <- names.row <- NULL
    counter <- 0
    message <- ""
    num.pieces <- length(unlist(nsnps.list))
    for(chr in names(nsnps.list)){
      k <- length(nsnps.list[[chr]])
      for(piece in 1:k){
        
        ## reference sample ##
        
        LD_rda_file <- LD.files[grep(x = LD.files, pattern = paste0("chr",chr,".",piece, "[\\._].*rda"))]
        LD_bim_file <- LD.files[grep(x = LD.files, pattern = paste0("chr",chr,".",piece, "[\\._].*bim"))]
        load(file=paste(LD.path, LD_rda_file, sep = "/"))
        snps.ref.df <- read.table(paste(LD.path, LD_bim_file, sep = "/"))
        
        colnames(snps.ref.df) <- c("chr","id","non","pos","A1","A2")
        snps.ref <- snps.ref.df$id
        A2.ref <- snps.ref.df$A2
        names(A2.ref) <- snps.ref
        
        gwas1.df.subset <- gwas1.df %>% filter(SNP %in% snps.ref) %>% distinct(SNP, A1, A2, .keep_all = TRUE)
        
        ## Check if there are multiallelic or duplicated SNPs
        if(any(duplicated(gwas1.df.subset$SNP)) == TRUE){
          gwas1.df.subset.duplicated <- gwas1.df.subset %>% 
            mutate(row.num = 1:n()) %>% 
            filter(SNP == SNP[duplicated(SNP)]) %>%
            mutate(SNP_A1_A2 = paste(SNP, A1, A2, sep = "_"))
          snps.ref.df.duplicated <- snps.ref.df %>%
            filter(id %in% gwas1.df.subset.duplicated$SNP)
          SNP_A1_A2.valid <- c(paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A1, snps.ref.df.duplicated$A2, sep = "_"),
                               paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A2, snps.ref.df.duplicated$A1, sep = "_"))
          row.remove <- gwas1.df.subset.duplicated %>% filter(!(SNP_A1_A2 %in% SNP_A1_A2.valid)) %>% select(row.num) %>% unlist()
          gwas1.df.subset <- gwas1.df.subset[-row.remove,]
        }
        
        bhat1.raw <- gwas1.df.subset[, "Z"] / sqrt(gwas1.df.subset[, "N"])
        A2.gwas1 <- gwas1.df.subset[, "A2"]
        names(bhat1.raw) <- names(A2.gwas1) <- gwas1.df.subset$SNP
        idx.sign1 <- A2.gwas1 == A2.ref[names(A2.gwas1)]
        bhat1.raw <- bhat1.raw*(2*as.numeric(idx.sign1)-1)
        
        
        
        gwas2.df.subset <- gwas2.df %>% filter(SNP %in% snps.ref) %>% distinct(SNP, A1, A2, .keep_all = TRUE)
        
        ## Check if there are multiallelic or duplicated SNPs
        if(any(duplicated(gwas2.df.subset$SNP)) == TRUE){
          gwas2.df.subset.duplicated <- gwas2.df.subset %>% 
            mutate(row.num = 1:n()) %>% 
            filter(SNP == SNP[duplicated(SNP)]) %>%
            mutate(SNP_A1_A2 = paste(SNP, A1, A2, sep = "_"))
          snps.ref.df.duplicated <- snps.ref.df %>%
            filter(id %in% gwas2.df.subset.duplicated$SNP)
          SNP_A1_A2.valid <- c(paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A1, snps.ref.df.duplicated$A2, sep = "_"),
                               paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A2, snps.ref.df.duplicated$A1, sep = "_"))
          row.remove <- gwas2.df.subset.duplicated %>% filter(!(SNP_A1_A2 %in% SNP_A1_A2.valid)) %>% select(row.num) %>% unlist()
          gwas2.df.subset <- gwas2.df.subset[-row.remove,]
        }
        bhat2.raw <- gwas2.df.subset[, "Z"] / sqrt(gwas2.df.subset[, "N"])
        A2.gwas2 <- gwas2.df.subset[, "A2"]
        names(bhat2.raw) <- names(A2.gwas2) <- gwas2.df.subset$SNP
        idx.sign2 <- A2.gwas2 == A2.ref[names(A2.gwas2)]
        bhat2.raw <- bhat2.raw*(2*as.numeric(idx.sign2)-1)
        
        M <- length(LDsc)
        bhat1 <- bhat2 <- numeric(M)
        names(bhat1) <- names(bhat2) <- snps.ref
        bhat1[names(bhat1.raw)] <- bhat1.raw
        bhat2[names(bhat2.raw)] <- bhat2.raw
        
        
        a11 <- bhat1**2
        a22 <- bhat2**2
        a12 <- bhat1*bhat2
        
        reg = lm(a11~ LDsc)
        h11.ols <- c(summary(reg)$coef[1:2,1:2]*c(N1,M))
        
        reg = lm(a22~ LDsc)
        h22.ols <- c(summary(reg)$coef[1:2,1:2]*c(N2,M))
        
        reg = lm(a12~ LDsc)
        if (N0>0) h12.ols = c(summary(reg)$coef[1:2,1:2]*c((N0/p1/p2),M))
        if (N0==0) h12.ols = c(summary(reg)$coef[1:2,1:2]*c(N,M))
        
        ##  ................................ weighted LS: use estimated h2
        ## vars from Bulik-Sullivan
        
        h11v = (h11.ols[2]*LDsc/M + 1/N1)^2
        h22v = (h22.ols[2]*LDsc/M + 1/N2)^2
        
        reg = lm(a11~ LDsc, weight=1/h11v)
        h11.wls <- c(summary(reg)$coef[1:2,1:2]*c(N1,M))
        
        reg = lm(a22~ LDsc, weight=1/h22v)
        h22.wls <- c(summary(reg)$coef[1:2,1:2]*c(N2,M))
        
        if (N0>0) h12v = sqrt(h11v*h22v) + (h12.ols[2]*LDsc/M + p1*p2*rho12/N0)^2
        if (N0==0) h12v = sqrt(h11v*h22v) + (h12.ols[2]*LDsc/M)^2
        
        reg = lm(a12~ LDsc, weight=1/h12v)
        if (N0>0) h12.wls = c(summary(reg)$coef[1:2,1:2]*c((N0/p1/p2),M))
        if (N0==0) h12.wls = c(summary(reg)$coef[1:2,1:2]*c(N,M))
        
        ## .................................  likelihood based
        ## ....  estimate h2s
        bstar1 = crossprod(V,bhat1)  ## 
        bstar2 = crossprod(V,bhat2)  ##
        
        opt = optim(c(h11.wls[2],1), llfun, N=N1, Nref=Nref, lam=lam, bstar=bstar1, M=M,
                    lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        
        h11.hdl = opt$par
        
        opt = optim(c(h22.wls[2],1), llfun, N=N2, Nref=Nref, lam=lam, bstar=bstar2, M=M,
                    lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h22.hdl = opt$par
        
        opt=  optim(c(h12.wls[2],rho12), llfun.gcov.part.2, h11=h11.hdl, h22=h22.hdl, 
                    rho12=rho12, M=M, N1=N1, N2=N2, N0=N0, Nref=Nref, 
                    lam0=lam, lam1=lam, lam2=lam, 
                    bstar1=bstar1, bstar2=bstar2,
                    lim=lim, method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
        h12.hdl = opt$par
        
        HDL11.df <- rbind(HDL11.df, h11.hdl)
        HDL22.df <- rbind(HDL22.df, h22.hdl)
        HDL12.df <- rbind(HDL12.df, h12.hdl)
        
        bstar1.v <- c(bstar1.v, list(bstar1))
        bstar2.v <- c(bstar2.v, list(bstar2))
        lam.v <- c(lam.v, list(lam))
        
        ## Report progress ##
        
        counter <- counter + 1
        value <- round(counter/num.pieces*100)
        backspaces <- paste(rep("\b", nchar(message)), collapse = "")
        message <- paste("Estimation is ongoing ... ", value, "%", sep = "", 
                         collapse = "")
        cat(backspaces, message, sep = "")
        
      }
    }
    cat("\n")
    rownames(HDL11.df) <- rownames(HDL22.df) <- rownames(HDL12.df) <- names.row
    
    h1_2 <- sum(HDL11.df[,1])
    h2_2 <- sum(HDL22.df[,1])
    gen.cov <- sum(HDL12.df[,1])
    
    
    ##### Estimated likelihood for h12 + h11,h22 independent #####
    
    cat("\n")
    cat("Integrating piecewise results \n")
    if(output.file != ""){
      cat("\n", file = output.file, append = T)
      cat("Integrating piecewise results \n", file = output.file, append = T)
    }
    M.ref <- sum(unlist(nsnps.list))
    eigen_select.fun <- function(x,k){
      return(x[1:k])
    }
    
    eigen_select_num.fun <- function(eigen.percent, cut.percent){
      if(!(eigen.percent[length(eigen.percent)] > cut.percent)){
        return(length(eigen.percent))
      } else{
        return(which(eigen.percent > cut.percent)[1])
      }
    }
    
    
    
    if(eigen.cut == "automatic"){
      eigen.num.v.90 <- eigen.num.v.95 <- eigen.num.v.99 <- c()
      nsnps.v <- unlist(nsnps.list)
      for(i in 1:length(nsnps.v)){
        lam.i <- lam.v[[i]]
        eigen.percent <- numeric(length(lam.i))
        temp <- 0
        for(j in 1:length(lam.i)){
          temp <- temp + lam.i[j]
          eigen.percent[j] <- temp/nsnps.v[i]
        }
        eigen.num.90 <- eigen_select_num.fun(eigen.percent, 0.9)
        eigen.num.95 <- eigen_select_num.fun(eigen.percent, 0.95)
        eigen.num.99 <- eigen_select_num.fun(eigen.percent, 0.99)
        
        eigen.num.v.90 <- c(eigen.num.v.90, eigen.num.90)
        eigen.num.v.95 <- c(eigen.num.v.95, eigen.num.95)
        eigen.num.v.99 <- c(eigen.num.v.99, eigen.num.99)
      }
      
      lam.v.90 <- mapply(eigen_select.fun,lam.v,eigen.num.v.90)
      bstar1.v.90 <- mapply(eigen_select.fun,bstar1.v,eigen.num.v.90)
      bstar2.v.90 <- mapply(eigen_select.fun,bstar2.v,eigen.num.v.90)
      
      opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.90), bstar=unlist(bstar1.v.90), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h11.hdl.90 = opt$par
      opt = optim(c(h2_2,1), llfun, N=N2, Nref=Nref, lam=unlist(lam.v.90), bstar=unlist(bstar2.v.90), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h22.hdl.90 = opt$par
      
      if(sum(unlist(eigen.num.v.90)) == sum(unlist(eigen.num.v.95))){
        lam.v.use <- lam.v.90
        bstar1.v.use <- bstar1.v.90
        bstar2.v.use <- bstar2.v.90
        h11.hdl.use <- h11.hdl.90
        h22.hdl.use <- h22.hdl.90
        eigen.use <- 0.9
      } else{
        lam.v.95 <- mapply(eigen_select.fun,lam.v,eigen.num.v.95)
        bstar1.v.95 <- mapply(eigen_select.fun,bstar1.v,eigen.num.v.95)
        bstar2.v.95 <- mapply(eigen_select.fun,bstar2.v,eigen.num.v.95)
        
        opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.95), bstar=unlist(bstar1.v.95), M=M.ref,
                    lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h11.hdl.95 = opt$par
        opt = optim(c(h2_2,1), llfun, N=N2, Nref=Nref, lam=unlist(lam.v.95), bstar=unlist(bstar2.v.95), M=M.ref,
                    lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h22.hdl.95 = opt$par
        
        if(sum(unlist(eigen.num.v.95)) == sum(unlist(eigen.num.v.99))){
          if(h11.hdl.95[1] != 0 &&
             h22.hdl.95[1] != 0 &&
             (h11.hdl.90[1]-h11.hdl.95[1])/abs(h11.hdl.95[1]) < 0.2 &&
             (h22.hdl.90[1]-h22.hdl.95[1])/abs(h22.hdl.95[1]) < 0.2){
            lam.v.use <- lam.v.95
            bstar1.v.use <- bstar1.v.95
            bstar2.v.use <- bstar2.v.95
            h11.hdl.use <- h11.hdl.95
            h22.hdl.use <- h22.hdl.95
            eigen.use <- 0.95
          } else{
            lam.v.use <- lam.v.90
            bstar1.v.use <- bstar1.v.90
            bstar2.v.use <- bstar2.v.90
            h11.hdl.use <- h11.hdl.90
            h22.hdl.use <- h22.hdl.90
            eigen.use <- 0.9
          }
        } else{
          lam.v.99 <- mapply(eigen_select.fun,lam.v,eigen.num.v.99)
          bstar1.v.99 <- mapply(eigen_select.fun,bstar1.v,eigen.num.v.99)
          bstar2.v.99 <- mapply(eigen_select.fun,bstar2.v,eigen.num.v.99)
          
          opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.99), bstar=unlist(bstar1.v.99), M=M.ref,
                      lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
          h11.hdl.99 = opt$par
          opt = optim(c(h2_2,1), llfun, N=N2, Nref=Nref, lam=unlist(lam.v.99), bstar=unlist(bstar2.v.99), M=M.ref,
                      lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
          h22.hdl.99 = opt$par
          
          if(h11.hdl.99[1] != 0 &&
             h22.hdl.99[1] != 0 &&
             (h11.hdl.90[1]-h11.hdl.99[1])/abs(h11.hdl.99[1]) < 0.2 &&
             (h22.hdl.90[1]-h22.hdl.99[1])/abs(h22.hdl.99[1]) < 0.2){
            lam.v.use <- lam.v.99
            bstar1.v.use <- bstar1.v.99
            bstar2.v.use <- bstar2.v.99
            h11.hdl.use <- h11.hdl.99
            h22.hdl.use <- h22.hdl.99
            eigen.use <- 0.99
          } else{
            if(h11.hdl.95[1] != 0 &&
               h22.hdl.95[1] != 0 &&
               (h11.hdl.90[1]-h11.hdl.95[1])/abs(h11.hdl.95[1]) < 0.2 &&
               (h22.hdl.90[1]-h22.hdl.95[1])/abs(h22.hdl.95[1]) < 0.2){
              lam.v.use <- lam.v.95
              bstar1.v.use <- bstar1.v.95
              bstar2.v.use <- bstar2.v.95
              h11.hdl.use <- h11.hdl.95
              h22.hdl.use <- h22.hdl.95
              eigen.use <- 0.95
            } else{
              lam.v.use <- lam.v.90
              bstar1.v.use <- bstar1.v.90
              bstar2.v.use <- bstar2.v.90
              h11.hdl.use <- h11.hdl.90
              h22.hdl.use <- h22.hdl.90
              eigen.use <- 0.9
            }
          }
        }
      }
    } else if(!is.na(as.numeric(eigen.cut))){
      eigen.cut <- as.numeric(eigen.cut)
      eigen.use <- eigen.cut
      eigen.num.v.cut <- c()
      nsnps.v <- unlist(nsnps.list)
      for(i in 1:length(nsnps.v)){
        lam.i <- lam.v[[i]]
        eigen.percent <- numeric(length(lam.i))
        temp <- 0
        for(j in 1:length(lam.i)){
          temp <- temp + lam.i[j]
          eigen.percent[j] <- temp/nsnps.v[i]
        }
        eigen.num.cut <- eigen_select_num.fun(eigen.percent, eigen.cut)
        eigen.num.v.cut <- c(eigen.num.v.cut, eigen.num.cut)
      }
      
      lam.v.cut <- mapply(eigen_select.fun,lam.v,eigen.num.v.cut)
      bstar1.v.cut <- mapply(eigen_select.fun,bstar1.v,eigen.num.v.cut)
      bstar2.v.cut <- mapply(eigen_select.fun,bstar2.v,eigen.num.v.cut)
      
      opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.cut), bstar=unlist(bstar1.v.cut), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h11.hdl.cut = opt$par
      opt = optim(c(h2_2,1), llfun, N=N2, Nref=Nref, lam=unlist(lam.v.cut), bstar=unlist(bstar2.v.cut), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h22.hdl.cut = opt$par
      
      lam.v.use <- lam.v.cut
      bstar1.v.use <- bstar1.v.cut
      bstar2.v.use <- bstar2.v.cut
      
      h11.hdl.use <- h11.hdl.cut
      h22.hdl.use <- h22.hdl.cut
    } 
    
    h11 <- h11.hdl.use[1]
    h22 <- h22.hdl.use[1]
    opt=  optim(c(gen.cov,rho12), llfun.gcov.part.2, h11=h11.hdl.use, h22=h22.hdl.use,
                rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                lam0=unlist(lam.v.use), lam1=unlist(lam.v.use), lam2=unlist(lam.v.use),
                bstar1=unlist(bstar1.v.use), bstar2=unlist(bstar2.v.use),
                lim=lim, method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
    if(opt$convergence != 0){
      starting.value.v <- c(0,-sqrt(h11*h22)*0.5, sqrt(h11*h22)*0.5)
      k <- 1
      while(opt$convergence != 0){
        starting.value <- starting.value.v[k]
        opt=  optim(c(starting.value,rho12), llfun.gcov.part.2, h11=h11.hdl.use, h22=h22.hdl.use,
                    rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                    lam0=unlist(lam.v.use), lam1=unlist(lam.v.use), lam2=unlist(lam.v.use),
                    bstar1=unlist(bstar1.v.use), bstar2=unlist(bstar2.v.use),
                    lim=lim, method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
        k <- k + 1
        if(k > length(starting.value.v)){
          error.message <- "Algorithm failed to converge after trying different initial values. \n"
          if(output.file != ""){
            cat(error.message, file = output.file, append = T)
          }
          stop(error.message)
        }
      }
    }
    h12.hdl.use <- opt$par
    h12 <- h12.hdl.use[1]
    rg <- h12.hdl.use[1]/sqrt(h11.hdl.use[1]*h22.hdl.use[1])
    
    if(intercept.output == T){
      h11.intercept <- h11.hdl.use[2]
      h22.intercept <- h22.hdl.use[2]
      h12.intercept <- h12.hdl.use[2]
    }
    
    output <- function(value){
      if(is.na(value)){
        value.out <- NA
      } else if(abs(value) < 1e-4){
        value.out <- formatC(value, format = "e", digits = 2)
      } else {
        value.out <- round(value, digits = 4)
      }
    }
    
    cat("\n")
    cat("Point estimates: \n")
    cat("Heritability of phenotype 1: ", 
        output(h11), 
        "\n")
    cat("Heritability of phenotype 2: ", 
        output(h22), 
        "\n")
    cat("Genetic Covariance: ", 
        output(h12), 
        "\n")
    cat("Genetic Correlation: ", 
        output(rg), 
        "\n")
    if(h11 == 0 | h22 == 0 ){
      cat("Warning: Heritability of one trait was estimated to be 0, which may due to:
          1) The true heritability is very small;
          2) The sample size is too small;
          3) Many SNPs in the chosen reference panel misses in the GWAS;
          4) There is a severe mismatch between the GWAS population and the population for computing reference panel")
    }
    cat("\n")
    
    
    cat("Continuing computing standard error with jackknife \n")
    if(output.file != ""){
      cat("Continuing computing standard error with jackknife \n", file = output.file, append = T)
    }
    counter <- 0
    message <- ""
    rg.jackknife <- h11.jackknife <- h12.jackknife <- h22.jackknife <- numeric(length(lam.v))
    if(intercept.output == T){
      h11.intercept.jackknife <- h12.intercept.jackknife <- h22.intercept.jackknife <- numeric(length(lam.v))
    }
    for(i in 1:length(lam.v)){
      opt = optim(h11.hdl.use, llfun, N=N1, Nref=Nref, lam=unlist(lam.v.use[-i]), bstar=unlist(bstar1.v.use[-i]), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h11.hdl.jackknife = opt$par
      
      
      opt = optim(h22.hdl.use, llfun, N=N2, Nref=Nref, lam=unlist(lam.v.use[-i]), bstar=unlist(bstar2.v.use[-i]), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h22.hdl.jackknife = opt$par
      
      opt=  optim(h12.hdl.use, llfun.gcov.part.2, h11=h11.hdl.use, h22=h22.hdl.use,
                  rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                  lam0=unlist(lam.v.use[-i]), lam1=unlist(lam.v.use[-i]), lam2=unlist(lam.v.use[-i]),
                  bstar1=unlist(bstar1.v.use[-i]), bstar2=unlist(bstar2.v.use[-i]),
                  lim=lim, method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
      h12.hdl.jackknife = opt$par
      h11.jackknife[i] <- h11.hdl.jackknife[1]
      h12.jackknife[i] <- h12.hdl.jackknife[1]
      h22.jackknife[i] <- h22.hdl.jackknife[1]
      rg.jackknife[i] <- h12.hdl.jackknife[1]/sqrt(h11.hdl.jackknife[1]*h22.hdl.jackknife[1])
      
      if(intercept.output == T){
        h11.intercept.jackknife[i] <- h11.hdl.jackknife[2]
        h12.intercept.jackknife[i] <- h12.hdl.jackknife[2]
        h22.intercept.jackknife[i] <- h22.hdl.jackknife[2]
      }
      
      ## Report progress ##
      
      counter <- counter + 1
      value <- round(counter/length(lam.v)*100)
      backspaces <- paste(rep("\b", nchar(message)), collapse = "")
      message <- paste("Progress... ", value, "%", sep = "", 
                       collapse = "")
      cat(backspaces, message, sep = "")
    }
    rg.jackknife <- rg.jackknife[!is.infinite(rg.jackknife)]
    h11.se <-  sqrt(mean((h11.jackknife - mean(h11.jackknife))^2)*(length(h11.jackknife) - 1))
    h12.se <-  sqrt(mean((h12.jackknife - mean(h12.jackknife))^2)*(length(h12.jackknife) - 1))
    h22.se <-  sqrt(mean((h22.jackknife - mean(h22.jackknife))^2)*(length(h22.jackknife) - 1))
    rg.se <-  sqrt(mean((rg.jackknife - mean(rg.jackknife))^2)*(length(rg.jackknife) - 1))
    P <- pchisq((rg/rg.se)^2, df = 1, lower.tail = FALSE)
    
    if(intercept.output == FALSE){
      estimates.df <- matrix(c(h11,h22,h12,rg, 
                               h11.se,h22.se,h12.se,rg.se),
                             nrow=4,ncol=2)
      rownames(estimates.df) <- c("Heritability_1", "Heritability_2", "Genetic_Covariance", "Genetic_Correlation")
      colnames(estimates.df) <- c("Estimate", "se")
    } else{
      h11.intercept.se <-  sqrt(mean((h11.intercept.jackknife - mean(h11.intercept.jackknife))^2)*(length(h11.intercept.jackknife) - 1))
      h12.intercept.se <-  sqrt(mean((h12.intercept.jackknife - mean(h12.intercept.jackknife))^2)*(length(h12.intercept.jackknife) - 1))
      h22.intercept.se <-  sqrt(mean((h22.intercept.jackknife - mean(h22.intercept.jackknife))^2)*(length(h22.intercept.jackknife) - 1))
      estimates.df <- matrix(c(h11,h22,h12,rg, h11.intercept, h22.intercept, h12.intercept,
                               h11.se,h22.se,h12.se,rg.se, h11.intercept.se, h22.intercept.se, h12.intercept.se),
                             nrow=7,ncol=2)
      rownames(estimates.df) <- c("Heritability_1", "Heritability_2", "Genetic_Covariance", "Genetic_Correlation", 
                                  "Heritability_1_intercept", "Heritability_2_intercept", "Genetic_Covariance_intercept")
      colnames(estimates.df) <- c("Estimate", "se")
    }
    
    
    
    
    if(is.na(P)){
      p.out <- NA
    } else{
      p.out <- formatC(P, format = "e", digits = 2)
    }
    
    end.time <- date()
    cat("\n")
    cat("\n")
    cat("Heritability of phenotype 1: ", 
        output(h11), 
        paste0("(", output(h11.se), ") \n"))
    cat("Heritability of phenotype 2: ", 
        output(h22), 
        paste0("(", output(h22.se), ") \n"))
    cat("Genetic Covariance: ", 
        output(h12), 
        paste0("(", output(h12.se), ") \n"))
    cat("Genetic Correlation: ", 
        output(rg), 
        paste0("(", output(rg.se), ") \n"))
    cat("P: ",p.out,"\n")
    if(h11 == 0 | h22 == 0 ){
      cat("Warning: Heritability of one trait was estimated to be 0, which may due to:
          1) The true heritability is very small;
          2) The sample size is too small;
          3) Many SNPs in the chosen reference panel misses in the GWAS.")
    }
    cat("\n")
    cat("Analysis finished at",end.time,"\n")
    
    if(output.file != ""){
      cat("\n", file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
      cat("Heritability of phenotype 1: ", 
          output(h11), 
          paste0("(", output(h11.se), ") \n"), file = output.file, append = TRUE)
      cat("Heritability of phenotype 2: ", 
          output(h22), 
          paste0("(", output(h22.se), ") \n"), file = output.file, append = TRUE)
      cat("Genetic Covariance: ", 
          output(h12), 
          paste0("(", output(h12.se), ") \n"), file = output.file, append = TRUE)
      cat("Genetic Correlation: ", 
          output(rg), 
          paste0("(", output(rg.se), ") \n"), file = output.file, append = TRUE)
      cat("P: ",p.out,"\n", file = output.file, append = TRUE)
      if(h11 == 0 | h22 == 0 ){
        cat("Warning: Heritability of one trait was estimated to be 0, which may due to:
            1) The true heritability is very low;
            2) The sample size of the GWAS is too small;
            3) Many SNPs in the chosen reference panel are absent in the GWAS;
            4) There is a severe mismatch between the GWAS population and the population for computing reference panel", file = output.file, append = TRUE)
      }
      cat("\n", file = output.file, append = TRUE)
      cat("Analysis finished at",end.time,"\n", file = output.file, append = TRUE)
      cat("The results were saved to", output.file)
      cat("\n")
      cat("The results were saved to", output.file, file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
    }
    
    if(jackknife.df == TRUE){
      jackknife.df <- rbind(h11.jackknife, h22.jackknife, h12.jackknife, rg.jackknife)
      rownames(jackknife.df) <- c("Heritability_1", "Heritability_2", "Genetic_Covariance", "Genetic_Correlation")
      return(list(rg = rg, rg.se = rg.se, P = P, estimates.df = estimates.df, eigen.use = eigen.use, jackknife.df = jackknife.df))
    }
    
    return(list(rg = rg, rg.se = rg.se, P = P, estimates.df = estimates.df, eigen.use = eigen.use))
  }
