 ## Deal with mismatch between gwas snp and reference snp
        
        match.gwas1 <- match(snps.ref, gwas1.df$SNP)
        snps.miss.gwas1.idx <- which(is.na(match.gwas1))
        snps.miss.gwas1 <- snps.ref[snps.miss.gwas1.idx]
        snps.nonmiss.gwas1.idx <- which(!is.na(match.gwas1))
        gwas1.df.subset <- gwas1.df %>% slice(match.gwas1)
        
        match.gwas2 <- match(snps.ref, gwas2.df$SNP)
        snps.miss.gwas2.idx <- which(is.na(match.gwas2))  
        snps.miss.gwas2 <- snps.ref[snps.miss.gwas2.idx]
        snps.nonmiss.gwas2.idx <- which(!is.na(match.gwas2))
        gwas2.df.subset <- gwas2.df %>% slice(match.gwas2)
        

        ## standardized bhat
        bhat1 <- gwas1.df.subset[, "Z"] / sqrt(gwas1.df.subset[, "N"])
        bhat2 <- gwas2.df.subset[, "Z"] / sqrt(gwas2.df.subset[, "N"])
        A2.gwas1 <- gwas1.df.subset[, "A2"]
        A2.gwas2 <- gwas2.df.subset[, "A2"]
        
        
        M <- length(LDsc)
        bhat1.full <- bhat2.full <- numeric(M)

        LDsc1 <- LDsc[setdiff(1:length(LDsc),snps.miss.gwas1.idx)]
        idx1 <- A2.gwas1 == A2.ref[snps.nonmiss.gwas1.idx]

        LDsc2 <- LDsc[setdiff(1:length(LDsc),snps.miss.gwas2.idx)]
        idx2 <- A2.gwas2 == A2.ref[snps.nonmiss.gwas2.idx]

        bhat1 <- bhat1*(2*as.numeric(idx1)-1)
        bhat2 <- bhat2*(2*as.numeric(idx2)-1)
        bhat1.full[snps.nonmiss.gwas1.idx] <- bhat1
        bhat2.full[snps.nonmiss.gwas2.idx] <- bhat2
        
        a11 <- bhat1**2
        a22 <- bhat2**2
        
        snp.gwas1.subset.both.idx <- setdiff(1:length(bhat1),match(snps.miss.gwas2, gwas1.df.subset$SNP))  # index of snps in gwas1.df.subset, which is also in gwas2
        snp.gwas2.subset.both.idx <- setdiff(1:length(bhat2),match(snps.miss.gwas1, gwas2.df.subset$SNP))  # index of snps in gwas1.df.subset, which is also in gwas2
        a12 <- bhat1[snp.gwas1.subset.both.idx]*bhat2[snp.gwas2.subset.both.idx]
        LDsc12 <- LDsc[setdiff(1:length(LDsc),c(snps.miss.gwas1.idx, snps.miss.gwas2.idx))]
        
        reg = lm(a11~ LDsc1)
        h11.ols <- c(summary(reg)$coef[1:2,1:2]*c(N1,M))
        
        reg = lm(a22~ LDsc2)
        h22.ols <- c(summary(reg)$coef[1:2,1:2]*c(N2,M))
        
        reg = lm(a12~ LDsc12)
        if (N0>0) h12.ols = c(summary(reg)$coef[1:2,1:2]*c((N0/p1/p2),M))
        if (N0==0) h12.ols = c(summary(reg)$coef[1:2,1:2]*c(N,M))
        
        ##  ................................ weighted LS: use estimated h2
        ## vars from Bulik-Sullivan
        
        h11v = (h11.ols[2]*LDsc1/M + 1/N1)^2
        h22v = (h22.ols[2]*LDsc2/M + 1/N2)^2
        
        reg11 = lm(a11~ LDsc1, weight=1/h11v)
        h11.wls <- c(summary(reg11)$coef[1:2,1:2]*c(N1,M))
        
        reg22 = lm(a22~ LDsc2, weight=1/h22v)
        h22.wls <- c(summary(reg22)$coef[1:2,1:2]*c(N2,M))
        
        if (N0>0) h12v = sqrt(h11v[snp.gwas1.subset.both.idx]*h22v[snp.gwas2.subset.both.idx]) + (h12.ols[2]*LDsc12/M + p1*p2*rho12/N0)^2
        if (N0==0) h12v = sqrt(h11v[snp.gwas1.subset.both.idx]*h22v[snp.gwas2.subset.both.idx]) + (h12.ols[2]*LDsc12/M)^2
        
        reg12 = lm(a12~ LDsc12, weight=1/h12v)
        if (N0>0) h12.wls = c(summary(reg12)$coef[1:2,1:2]*c((N0/p1/p2),M))
        if (N0==0) h12.wls = c(summary(reg12)$coef[1:2,1:2]*c(N,M))
        
