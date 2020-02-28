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
#' @param Nref Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 336000.
#' @param N0 Number of individuals included in both cohorts. However, the estimated genetic correlation is usually robust against misspecified N0. 
#' If not given, the default value is set to the minimum sample size across all SNPs in cohort 1 and cohort 2.
#' @param output.file Where the log and results should be written. If you do not specify a file, the log will be printed on the console.
#' 
#' @note Users can download the eigenvalues and eigenvectors of LD correlation matrices from https://www.dropbox.com/sh/ai2o21gxklhlvs3/AABPD7nKv3nXcvmoQQ3cGh9Qa?dl=0. 
#' These are the LD matrices and their eigen-decomposition from 336,000 genomic British UK Biobank individuals. Two sets of reference panel are provided: 
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
#' }
#' 
#' @author Zheng Ning
#' 
#' @references 
#' Ning Z, Pawitan Y, Shen X (2019). High-definition likelihood inference of genetic correlations 
#' across human complex traits. \emph{Submitted}.
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
  function(gwas1.df, gwas2.df, LD.path, Nref = 336000, N0 = min(gwas1.df$N, gwas2.df$N), output.file = ""){
    
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
    
    if(file.exists(paste0(LD.path, "/overlap.snp.MAF.05.list.rda"))){
      load(file=paste0(LD.path, "/UKB_snp_counter_overlap_MAF_5.RData"))
      load(file=paste0(LD.path, "/overlap.snp.MAF.05.list.rda"))
    } else if(file.exists(paste0(LD.path, "UKB_snp_list_imputed.vector_form.RData"))){
      load(file=paste0(LD.path, "/UKB_snp_counter_imputed.RData"))
      load(file=paste0(LD.path, "/UKB_snp_list_imputed.vector_form.RData"))
      overlap.snp.MAF.05.list <- snps.list.imputed.vector
      nsnps.list <- nsnps.list.imputed
    } else{
      error.message <- "It seems this directory does not contain all files needed for HDL. Please check your LD.path again. The current version of HDL only support pre-computed LD reference panels."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
    
    k1 <- sum(gwas1.df$SNP %in% overlap.snp.MAF.05.list)
    k2 <- sum(gwas2.df$SNP %in% overlap.snp.MAF.05.list)
    k1.percent <- paste("(",round(100*k1 / length(overlap.snp.MAF.05.list), 2), "%)", sep="") 
    k2.percent <- paste("(",round(100*k2 / length(overlap.snp.MAF.05.list), 2), "%)", sep="") 

    cat(k1, "out of", length(overlap.snp.MAF.05.list), k1.percent, "SNPs in reference panel are available in GWAS 1."," \n")
    cat(k2, "out of", length(overlap.snp.MAF.05.list), k2.percent, "SNPs in reference panel are available in GWAS 2."," \n")
    if(output.file != ""){
      cat(k1, "out of", length(overlap.snp.MAF.05.list), k1.percent, "SNPs in reference panel are available in GWAS 1."," \n", file = output.file, append = T)
      cat(k2, "out of", length(overlap.snp.MAF.05.list), k2.percent, "SNPs in reference panel are available in GWAS 2."," \n", file = output.file, append = T)
    }
    if(k1 < length(overlap.snp.MAF.05.list)*0.99){
      error.message <- "More than 1% SNPs in reference panel are missed in GWAS 1. Please check."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
    if(k2 < length(overlap.snp.MAF.05.list)*0.99){
      error.message <- "More than 1% SNPs in reference panel are missed in GWAS 2. Please check."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
    
    ## stats
    N1 <- max(gwas1.df[, "N"])
    N2 <- max(gwas2.df[, "N"])
    N <- sqrt(N1)*sqrt(N2)
    
    ## samples for phenotypes
    p1 <- N0/N1
    p2 <- N0/N2
    rho12 <- cor(gwas1.df$Z, gwas2.df$Z, use = "complete.obs") 
    
    bstar1.v <- bstar2.v <- lam.v <- list()
    HDL11.df <- HDL12.df <- HDL22.df <- names.row <- NULL
    counter <- 0
    message <- ""
    num.pieces <- length(unlist(nsnps.list))
    for(chr in 1:22){
      k <- length(nsnps.list[[chr]])
      for(piece in 1:k){
        
        ## reference sample ##
        
        if(file.exists(paste0(LD.path, "/overlap.snp.MAF.05.list.rda"))){
          load(file=paste0(LD.path, "/ukb_chr",chr,".",piece,"_n336000_500banded_90eigen.rda"))
          snps.ref.df <- read.table(paste0(LD.path, "/ukb_chr",chr,".",piece,"_n336000.bim"))
        } else if(file.exists(paste0(LD.path, "/UKB_snp_list_imputed.vector_form.RData"))){
          load(file=paste0(LD.path, "/ukb_imputed_chr",chr,".",piece,"_n336000_500banded_99eigen.rda"))
          snps.ref.df <- read.table(paste0(LD.path, "/ukb_chr",chr,".",piece,"_n336000.imputed_clean.bim"))
        }
        colnames(snps.ref.df) <- c("chr","id","non","pos","A1","A2")
        snps.ref <- snps.ref.df$id
        A2.ref <- snps.ref.df$A2
        
        gwas1.df.subset <- gwas1.df %>% filter(SNP %in% snps.ref)
        gwas2.df.subset <- gwas2.df %>% filter(SNP %in% snps.ref)
        
        nonmissing1 <- which(snps.ref %in% gwas1.df.subset$SNP)
        nonmissing2 <- which(snps.ref %in% gwas2.df.subset$SNP)
        
        # Nonmissing SNPs 
        bhat1.0 <- gwas1.df.subset[, "Z"] / sqrt(gwas1.df.subset[, "N"])  
        bhat2.0 <- gwas2.df.subset[, "Z"] / sqrt(gwas2.df.subset[, "N"])  
        A2.gwas1 <- gwas1.df.subset[, "A2"]
        A2.gwas2 <- gwas2.df.subset[, "A2"]
        
        idx.sign1 <- A2.gwas1 == A2.ref[nonmissing1]
        idx.sign2 <- A2.gwas2 == A2.ref[nonmissing2]
        bhat1.0 <- bhat1.0*(2*as.numeric(idx.sign1)-1)
        bhat2.0 <- bhat2.0*(2*as.numeric(idx.sign2)-1)
        
        M <- length(LDsc)
        bhat1 <- bhat2 <- numeric(M)
        bhat1[nonmissing1] <- bhat1.0
        bhat2[nonmissing2] <- bhat2.0
        
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
                    lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        
        h11.hdl = opt$par
        
        opt = optim(c(h22.wls[2],1), llfun, N=N2, Nref=Nref, lam=lam, bstar=bstar2, M=M,
                    lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h22.hdl = opt$par
        
        opt=  optim(c(h12.wls[2],rho12), llfun.gcov.part.2, h11=h11.hdl, h22=h22.hdl, 
                    rho12=rho12, M=M, N1=N1, N2=N2, N0=N0, Nref=Nref, 
                    lam0=lam, lam1=lam, lam2=lam, 
                    bstar1=bstar1, bstar2=bstar2,
                    lim=exp(-18), method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
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
    opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v), bstar=unlist(bstar1.v), M=M.ref,
                lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
    h11.hdl = opt$par
    
    
    opt = optim(c(h2_2,1), llfun, N=N2, Nref=Nref, lam=unlist(lam.v), bstar=unlist(bstar2.v), M=M.ref,
                lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
    h22.hdl = opt$par
    
    opt=  optim(c(gen.cov,rho12), llfun.gcov.part.2, h11=h11.hdl, h22=h22.hdl,
                rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                lam0=unlist(lam.v), lam1=unlist(lam.v), lam2=unlist(lam.v),
                bstar1=unlist(bstar1.v), bstar2=unlist(bstar2.v),
                lim=exp(-18), method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
    h12.hdl = opt$par
    rg <- h12.hdl[1]/sqrt(h11.hdl[1]*h22.hdl[1])
    
    
    ##### Check whether rg is sensible #####
    counter <- 0
    message <- ""
    
    ## rg is sensible
    if(abs(rg) <1){
      cat("The estimated rg is within (-1,1). Continuing computing standard error. \n")
      if(output.file != ""){
        cat("The estimated rg is within (-1,1). Continuing computing standard error. \n", file = output.file, append = T)
      }
      set.seed(510)
      rg.jackknife <- length(lam.v)
      for(i in 1:length(lam.v)){
        opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v[-i]), bstar=unlist(bstar1.v[-i]), M=M.ref,
                    lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h11.hdl.jackknife = opt$par
        
        
        opt = optim(c(h2_2,1), llfun, N=N2, Nref=Nref, lam=unlist(lam.v[-i]), bstar=unlist(bstar2.v[-i]), M=M.ref,
                    lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h22.hdl.jackknife = opt$par
        
        opt=  optim(c(gen.cov,rho12), llfun.gcov.part.2, h11=h11.hdl, h22=h22.hdl,
                    rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                    lam0=unlist(lam.v[-i]), lam1=unlist(lam.v[-i]), lam2=unlist(lam.v[-i]),
                    bstar1=unlist(bstar1.v[-i]), bstar2=unlist(bstar2.v[-i]),
                    lim=exp(-18), method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
        h12.hdl.jackknife = opt$par
        rg.jackknife[i] <- h12.hdl.jackknife[1]/sqrt(h11.hdl.jackknife[1]*h22.hdl.jackknife[1])
        
        ## Report progress ##
        
        counter <- counter + 1
        value <- round(counter/length(lam.v)*100)
        backspaces <- paste(rep("\b", nchar(message)), collapse = "")
        message <- paste("Progress... ", value, "%", sep = "", 
                         collapse = "")
        cat(backspaces, message, sep = "")
        
      }
      rg.se <-  sqrt(mean((rg.jackknife - mean(rg.jackknife))^2)*(length(rg.jackknife) - 1))
      P <- pchisq((rg/rg.se)^2, df = 1, lower.tail = FALSE)
    } else{
      ## rg is not sensible, switch to full likelihood for rg estimation, slower
      cat("The estimated rg beyonds (-1,1). Switching to full likelihood. \n")
      if(output.file != ""){
        cat("The estimated rg beyonds (-1,1). Switching to full likelihood. \n", file = output.file, append = T)
      }
      opt=  optim(c(0.5,1, 0.5,1, sign(rg)*0.2, rho12), llfun.rg.full.likelihood, 
                  M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref, 
                  lam1=unlist(lam.v), lam2=unlist(lam.v),
                  bstar1=unlist(bstar1.v), bstar2=unlist(bstar2.v),
                  lim=exp(-18), method ='L-BFGS-B', lower=c(rep(0,4),-1,-10), upper=c(1, 10, 1, 10, 1, 10))
      if(opt$value > 1e5 | opt$convergence != 0){
        opt=  optim(c(0.5,1, 0.5,1, sign(rg)*0.5, rho12), llfun.rg.full.likelihood, 
                    M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref, 
                    lam1=unlist(lam.v), lam2=unlist(lam.v),
                    bstar1=unlist(bstar1.v), bstar2=unlist(bstar2.v),
                    lim=exp(-18), method ='L-BFGS-B', lower=c(rep(0,4),-1,-10), upper=c(1, 10, 1, 10, 1, 10))
      }
      if(opt$value > 1e5 | opt$convergence != 0){
        opt=  optim(c(0.5,1, 0.5,1, sign(rg)*0.8, rho12), llfun.rg.full.likelihood, 
                    M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref, 
                    lam1=unlist(lam.v), lam2=unlist(lam.v),
                    bstar1=unlist(bstar1.v), bstar2=unlist(bstar2.v),
                    lim=exp(-18), method ='L-BFGS-B', lower=c(rep(0,4),-1,-10), upper=c(1, 10, 1, 10, 1, 10))
      }
      if(opt$value > 1e5 | opt$convergence != 0){
        stop("Algorithm fails to converge. It may because heritabilities are too small or genetic correlation is very close to 1.")
      }
      rg.full.likelihood <- opt$par
      rg <- rg.full.likelihood[5]
      
      cat("rg is within (-1,1) now. Continuing computing standard error. \n")
      if(output.file != ""){
        cat("rg is within (-1,1) now. Continuing computing standard error. \n", file = output.file, append = T)
      }
      set.seed(510)
      rg.jackknife <- length(lam.v)
      for(i in 1:length(lam.v)){
        opt=  optim(rg.full.likelihood, llfun.rg.full.likelihood,
                    M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                    lam1=unlist(lam.v[-i]), lam2=unlist(lam.v[-i]),
                    bstar1=unlist(bstar1.v[-i]), bstar2=unlist(bstar2.v[-i]),
                    lim=exp(-18), method ='L-BFGS-B', lower=c(rep(0,4),-1,-10), upper=c(1, 10, 1, 10, 1, 10))
        rg.jackknife[i] <- opt$par[5]
        
        ## Report progress ##

        counter <- counter + 1
        value <- round(counter/length(lam.v)*100)
        backspaces <- paste(rep("\b", nchar(message)), collapse = "")
        message <- paste("Progress... ", value, "%", sep = "", 
                         collapse = "")
        cat(backspaces, message, sep = "")
        
      }
      rg.se <-  sqrt(mean((rg.jackknife - mean(rg.jackknife))^2)*(length(rg.jackknife) - 1))
      P <- pchisq((rg/rg.se)^2, df = 1, lower.tail = FALSE)
    }
    
    
      if(abs(rg) < 1e-4){
        rg.out <- formatC(rg, format = "e", digits = 2)
      } else{
        rg.out <- round(rg, digits = 4)
      }
      
      if(rg.se < 1e-4){
        rg.se.out <- formatC(rg.se, format = "e", digits = 2)
      } else{
        rg.se.out <- round(rg.se, digits = 4)
      }
      
      p.out <- formatC(P, format = "e", digits = 2)
      
      end.time <- date()
      cat("\n")
      cat("\n")
      cat("Genetic Correlation: ", 
          rg.out, 
          paste0("(", rg.se.out, ")"))
      cat("\n")
      cat("P: ",p.out)
      cat("\n")
      cat("\n")
      cat("Analysis finished at",end.time,"\n")
      
      if(output.file != ""){
      cat("\n", file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
      cat("Genetic Correlation: ", 
          rg.out, 
          paste0("(", rg.se.out, ")"), 
          file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
      cat("P: ",p.out, file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
      cat("Analysis finished at",end.time,"\n", file = output.file, append = TRUE)
      cat("The results were saved to", output.file)
      cat("\n")
      cat("The results were saved to", output.file, file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
    }
    
    return(list(rg = rg, rg.se = rg.se, P = P))
  }
