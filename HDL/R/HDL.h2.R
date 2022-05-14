#' High-definition likelihood inference of SNP-based heritabiity (HDL)
#' 
#' The function returns the estimate and standard error of the SNP-heritability of one trait based on GWAS summary statistics. 
#' 
#' @param gwas.df A data frame including GWAS summary statistics of genetic variants for a trait. 
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS. 
#' @param LD.path Path to the directory where linkage disequilibrium (LD) information is stored.
#' @param Nref Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 335265
#' @param output.file Where the log and results should be written. If you do not specify a file, the log will be printed on the console.
#' @param eigen.cut Which eigenvalues and eigenvectors in each LD score matrix should be used for HDL. 
#' Users are allowed to specify a numeric value between 0 and 1 for eigen.cut. For example, eigen.cut = 0.99 means using the leading eigenvalues explaining 99% of the variance
#' and their correspondent eigenvectors. If the default 'automatic' is used, the eigen.cut gives the most stable heritability estimates will be used. 
#' @param lim Tolerance limitation, default lim = exp(-18). 
#' @note Users can download the precomputed eigenvalues and eigenvectors of LD correlation matrices for European ancestry population. The download link can be found at https://github.com/zhenin/HDL/wiki/Reference-panels
#' These are the LD matrices and their eigen-decomposition from 335,265 genomic British UK Biobank individuals. Three sets of reference panel are provided: 
#' 1) 1,029,876 QCed UK Biobank imputed HapMap3 SNPs. The size is about 33 GB after unzipping. Although it takes more time, using the imputed panel provides more accurate estimates of genetic correlations. 
#' Therefore if the GWAS includes most of the HapMap3 SNPs, then we recommend using the imputed reference panel.
#' 2) 769,306 QCed UK Biobank imputed HapMap2 SNPs. The size is about 18 GB after unzipping.If one of your GWAS includes most of the HapMap 2 SNPs, but many SNPs (more than 1%) in the above HapMap 3 reference panel are absent, 
#' then this HapMap2 panel is more proper to be used for HDL. 
#' 3) 307,519 QCed UK Biobank Axiom Array SNPs. The size is about 7.5 GB after unzipping.
#' 
#' @return A list is returned with:
#' \itemize{
#' \item{h2 }{The estimated SNP-based heritability.}
#' \item{h2.se }{The standard error of the estimated heritability.}
#' \item{P }{P-value based on Wald test.}
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
#' ## The path to the directory where linkage disequilibrium (LD) information is stored.
#' LD.path <- "/Users/zhengning/Work/HDL/package/UKB_array_SVD_eigen90_extraction"
#' 
#' res.HDL <- HDL.h2(gwas1.example, LD.path)
#' res.HDL
#' }
#' @export
#' 

HDL.h2 <-
  function(gwas.df, LD.path, Nref = 335265, output.file = "", eigen.cut = "automatic", intercept.output = FALSE, fill.missing.N = NULL, lim = exp(-18)){
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
    
    gwas.df <- gwas.df %>% filter(SNP %in% snps.name.list)
    
    gwas.df$A1 <- toupper(as.character(gwas.df$A1))
    gwas.df$A2 <- toupper(as.character(gwas.df$A2))
    
    
    
    
    if(!("Z" %in% colnames(gwas.df))){
      if(("b" %in% colnames(gwas.df)) && ("se" %in% colnames(gwas.df))){
        if(abs(median(gwas.df$b) - 1) < 0.1){
          cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n")
          if(output.file != ""){
            cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n", file = output.file, append = T)
          }
          gwas.df$Z <- log(gwas.df$b) / gwas.df$se
        } else{
          gwas.df$Z <- gwas.df$b / gwas.df$se
        }
      } else{
        error.message <- "Z is not available in GWAS 1, meanwhile either b or se is missing. Please check."
        if(output.file != ""){
          cat(error.message, file = output.file, append = T)
        }
        stop(error.message)
        
      }
    }
    
    k1.0 <- length(unique(gwas.df$SNP))
    
    if(is.null(fill.missing.N)){
      gwas.df <- gwas.df %>% filter(!is.na(Z), !is.na(N))
      
    } else if(fill.missing.N == "min"){
      gwas.df <- gwas.df %>% filter(!is.na(Z))
      gwas.df$N[is.na(gwas.df$N)] <- min(gwas.df$N, na.rm = T)
    } else if(fill.missing.N == "max"){
      gwas.df <- gwas.df %>% filter(!is.na(Z))
      gwas.df$N[is.na(gwas.df$N)] <- max(gwas.df$N, na.rm = T)
    } else if(fill.missing.N == "median"){
      gwas.df <- gwas.df %>% filter(!is.na(Z))
      gwas.df$N[is.na(gwas.df$N)] <- median(gwas.df$N, na.rm = T)
    } else{
      error.message <- "If given, the argument fill.missing.N can only be one of below: 'min', 'max', 'median'."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
    
    k1 <- length(unique(gwas.df$SNP))
    k1.percent <- paste("(",round(100*k1 / length(snps.name.list), 2), "%)", sep="") 
    
    cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in the GWAS."," \n")
    if(output.file != ""){
      cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in the GWAS."," \n", file = output.file, append = T)
    }
    if(k1 < length(snps.name.list)*0.99){
      error.message <- "Warning: More than 1% SNPs in reference panel are missed in the GWAS. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      cat(error.message)
    }
    
    
    
    ## stats
    N1 <- median(gwas.df[, "N"])
    N <- N1
    
    
    bstar1.v <- lam.v <- list()
    HDL11.df <-names.row <- NULL
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
        
        gwas1.df.subset <- gwas.df %>% filter(SNP %in% snps.ref) %>% distinct(SNP, A1, A2, .keep_all = TRUE)
        
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
        
        
        
        M <- length(LDsc)
        bhat1 <- numeric(M)
        names(bhat1) <- snps.ref
        bhat1[names(bhat1.raw)] <- bhat1.raw
        
        a11 <- bhat1**2
        
        reg = lm(a11~ LDsc)
        h11.ols <- c(summary(reg)$coef[1:2,1:2]*c(N1,M))
        
        ##  ................................ weighted LS: use estimated h2
        ## vars from Bulik-Sullivan
        
        h11v = (h11.ols[2]*LDsc/M + 1/N1)^2
        
        reg = lm(a11~ LDsc, weight=1/h11v)
        h11.wls <- c(summary(reg)$coef[1:2,1:2]*c(N1,M))
        
        ## .................................  likelihood based
        ## ....  estimate h2s
        bstar1 = crossprod(V,bhat1)  ## 
        
        opt = optim(c(h11.wls[2],1), llfun, N=N1, Nref=Nref, lam=lam, bstar=bstar1, M=M,
                    lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        
        h11.hdl = opt$par
        
        
        HDL11.df <- rbind(HDL11.df, h11.hdl)
        
        bstar1.v <- c(bstar1.v, list(bstar1))
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
    rownames(HDL11.df) <- names.row
    h1_2 <- sum(HDL11.df[,1])
    
    
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
      
      opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.90), bstar=unlist(bstar1.v.90), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h11.hdl.90 = opt$par
      
      if(sum(unlist(eigen.num.v.90)) == sum(unlist(eigen.num.v.95))){
        lam.v.use <- lam.v.90
        bstar1.v.use <- bstar1.v.90
        h11.hdl.use <- h11.hdl.90
        eigen.use <- 0.9
      } else{
        lam.v.95 <- mapply(eigen_select.fun,lam.v,eigen.num.v.95)
        bstar1.v.95 <- mapply(eigen_select.fun,bstar1.v,eigen.num.v.95)
        
        opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.95), bstar=unlist(bstar1.v.95), M=M.ref,
                    lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h11.hdl.95 = opt$par
        
        if(sum(unlist(eigen.num.v.95)) == sum(unlist(eigen.num.v.99))){
          if(h11.hdl.95[1] != 0 &&
             (h11.hdl.90[1]-h11.hdl.95[1])/abs(h11.hdl.95[1]) < 0.2){
            lam.v.use <- lam.v.95
            bstar1.v.use <- bstar1.v.95
            h11.hdl.use <- h11.hdl.95
            eigen.use <- 0.95
          } else{
            lam.v.use <- lam.v.90
            bstar1.v.use <- bstar1.v.90
            h11.hdl.use <- h11.hdl.90
            eigen.use <- 0.9
          }
        } else{
          lam.v.99 <- mapply(eigen_select.fun,lam.v,eigen.num.v.99)
          bstar1.v.99 <- mapply(eigen_select.fun,bstar1.v,eigen.num.v.99)
          
          opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.99), bstar=unlist(bstar1.v.99), M=M.ref,
                      lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
          h11.hdl.99 = opt$par
          
          if(h11.hdl.99[1] != 0 &&
             (h11.hdl.90[1]-h11.hdl.99[1])/abs(h11.hdl.99[1]) < 0.2){
            lam.v.use <- lam.v.99
            bstar1.v.use <- bstar1.v.99
            h11.hdl.use <- h11.hdl.99
            eigen.use <- 0.99
          } else{
            if(h11.hdl.95[1] != 0 &&
               (h11.hdl.90[1]-h11.hdl.95[1])/abs(h11.hdl.95[1]) < 0.2){
              lam.v.use <- lam.v.95
              bstar1.v.use <- bstar1.v.95
              h11.hdl.use <- h11.hdl.95
              eigen.use <- 0.95
            } else{
              lam.v.use <- lam.v.90
              bstar1.v.use <- bstar1.v.90
              h11.hdl.use <- h11.hdl.90
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
      
      opt = optim(c(h1_2,1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.cut), bstar=unlist(bstar1.v.cut), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h11.hdl.cut = opt$par
      
      
      lam.v.use <- lam.v.cut
      bstar1.v.use <- bstar1.v.cut
      
      h11.hdl.use <- h11.hdl.cut
    } 
    
    h11 <- h11.hdl.use[1]
    
    if(intercept.output == T){
      h11.intercept <- h11.hdl.use[2]
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
    cat("Point estimate: \n")
    cat("Heritability: ", 
        output(h11), 
        "\n")
    if(h11 == 0){
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
    rg.jackknife <- h11.jackknife <- length(lam.v)
    if(intercept.output == T){
      h11.intercept.jackknife <- numeric(length(lam.v))
    }
    for(i in 1:length(lam.v)){
      opt = optim(h11.hdl.use, llfun, N=N1, Nref=Nref, lam=unlist(lam.v.use[-i]), bstar=unlist(bstar1.v.use[-i]), M=M.ref,
                  lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h11.hdl.jackknife = opt$par
      
      h11.jackknife[i] <- h11.hdl.jackknife[1]
      
      if(intercept.output == T){
        h11.intercept.jackknife[i] <- h11.hdl.jackknife[2]
      }
      
      ## Report progress ##
      
      counter <- counter + 1
      value <- round(counter/length(lam.v)*100)
      backspaces <- paste(rep("\b", nchar(message)), collapse = "")
      message <- paste("Progress... ", value, "%", sep = "", 
                       collapse = "")
      cat(backspaces, message, sep = "")
    }
    h11.se <-  sqrt(mean((h11.jackknife - mean(h11.jackknife))^2)*(length(h11.jackknife) - 1))
    
    if(intercept.output == TRUE){
      h11.intercept.se <-  sqrt(mean((h11.intercept.jackknife - mean(h11.intercept.jackknife))^2)*(length(h11.intercept.jackknife) - 1))
    }
    P <- pchisq((h11/h11.se)^2, df = 1, lower.tail = FALSE)
    
    if(is.na(P)){
      p.out <- NA
    } else{
      p.out <- formatC(P, format = "e", digits = 2)
    }
    
    end.time <- date()
    cat("\n")
    cat("\n")
    cat("Heritability: ", 
        output(h11), 
        paste0("(", output(h11.se), ") \n"))
    cat("P: ",p.out,"\n")
    if(h11 == 0){
      cat("Warning: Heritability of the trait was estimated to be 0, which may due to:
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
      cat("P: ",p.out,"\n", file = output.file, append = TRUE)
      if(h11 == 0){
        cat("Warning: Heritability of the trait was estimated to be 0, which may due to:
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
    if(intercept.output == TRUE){
      return(list(h2 = h11, h2.se = h11.se, intercept = h11.intercept, intercept.se = h11.intercept.se, P = P, eigen.use = eigen.use))
    }
    return(list(h2 = h11, h2.se = h11.se, P = P, eigen.use = eigen.use))
  }
