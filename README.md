What is HDL?
------------

High-Definition Likelihood (HDL) is a likelihood-based method for estimating genetic correlation using GWAS summary statistics. It reduces the variance of a genetic correlation estimate by about 50%, compared to LD Score regression (LDSC). Here, we provide a computational tool `HDL` to implement our method. Although `HDL` is written in R, you can use it with command line. So no worry if you are not an R user.

What data are required?
-----------------------

-   `gwas1.df` A data frame including GWAS summary statistics of genetic variants for trait 1. The input data frame should include following columns: `SNP`, SNP ID; `A1`, effect allele; `A2`, reference allele; `N`, sample size; `Z`, z-score; If `Z` is not given, alternatively, you may provide: `b`, estimate of marginal effect in GWAS; `se`, standard error of the estimates of marginal effects in GWAS. The summary statistics should look like this:

<!-- -->

    ##          SNP A1 A2      N        b       se      Z
    ## 1  rs3131962  G  A 205475 0.001004 0.004590 0.2187
    ## 2 rs12562034  A  G 205475 0.005382 0.005011 1.0740
    ## 3 rs11240779  A  G 205475 0.002259 0.003691 0.6119
    ## 4 rs57181708  G  A 205475 0.005401 0.005114 1.0562
    ## 5  rs4422948  G  A 205475 0.005368 0.003604 1.4893
    ## 6  rs4970383  A  C 205475 0.004685 0.003582 1.3080

-   `gwas2.df` A data frame including GWAS summary statistics of genetic variants for trait 2. The format is same as `gwas1.df`.

-   The eigenvalues and eigenvectors of LD matrices. For European population, we have computed the LD and their eigens from 336,000 ethnically British individuals in UK Biobank. You can download them [here](https://www.dropbox.com/sh/hkg8x55s03d9zbh/AACaGPDqpPjAkLzeanRjSLrNa?dl=0). The size is about 7.5 GB after unzipping. These eigens of LD are from eigen-decomposition results based on 307,519 QCed UK Biobank Axiom Array SNPs.

Installation
------------

Firstly you need an R (version &gt;= 3.5.2) in your computer. If you have not installed R, you can easily download and install it from <https://cran.r-project.org/>.

#### Command line user

In order to download the R package `HDL` and related scripts, you can clone `HDL` repository and then install it using `HDL.install.R`:

``` r
git clone https://github.com/zhenin/HDL
cd HDL
Rscript HDL.install.R
```

#### R user

You can use the following codes to install `HDL` from Github:

``` r
install.packages("https://github.com/zhenin/HDL/raw/master/HDL_1.0.tar.gz", repos=NULL)
```

and load the package via:

``` r
library(HDL)
```

Estimating genetic correlation using HDL
----------------------------------------

To illustrate how to use HDL, we include two cleaned UKB GWAS summary statistics datasets as examples. `gwas1.example.rda` is for birth weight; and `gwas2.example.rda` is for type 2 diabetes.

#### Command line user

Next, you can simply run `HDL.run.R` like below to use HDL:

``` r
Rscript HDL.run.R gwas1.df=gwas1.example.rda gwas2.df=gwas2.example.rda LD.path=UKB_SVD_eigen90_extraction output.file=test.Rout
```

There are several arguments you should pass to `HDL`. **Please note that when you specify arguments, there should not be any space on any side of `=`.**

-   Mandatory arguments
-   `gwas1.df`, path to the file of the GWAS results for trait 1;
-   `gwas2.df`, path to the file of the GWAS results for trait 2;
-   `LD.path`, path to the directory where linkage disequilibrium (LD) information is stored.
-   Optional arguments
-   `output.file`, where the log and results should be written. If you do not specify a file, the log will be printed in console;
-   `Nref`, sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 336000;
-   `N0`, number of individuals included in both cohorts. However, the estimated genetic correlation is usually robust against misspecified N0. If not given, the default value is set to the minimum sample size across all SNPs in cohort 1 and cohort 2.

#### R user

`HDL.rg` is the function to perform HDL. The arguments for `HDL.rg` is the same as above arguments for command line implementation. Let's have a try with example data as below

``` r
data(gwas1.example)
data(gwas2.example)
LD.path <- "/Users/zhengning/Work/HDL/package/UKB_SVD_eigen90_extraction"
res.HDL <- HDL.rg(gwas1.example, gwas2.example, LD.path)
res.HDL
```

A list is returned with

-   `rg`, the estimation of genetic correlation.
-   `rg.se`, the standard error of estimated genetic correlation.
-   `P`, the Wald test P-value for rg

For direct R documentation of `HDL.rg` function, you can use question mark in R:

``` r
?HDL.rg
```

Reading HDL results
-------------------

The first section provides version information of HDL package:

``` r
*********************************************************************
* High-definition likelihood (HDL)
* Version 1.0
*********************************************************************
```

The next section gives the proportion of overlap SNPs between GWAS summary statistics and reference panel. A low SNP overlap may lead to poor estimation.

``` r
Analysis starts 
307519 out of 307519 (100%) SNPs in reference panel are available in GWAS 1.  
307519 out of 307519 (100%) SNPs in reference panel are available in GWAS 2.  
```

The following section shows whether the estimated rg is within (-1,1). If not, we switch to an algorithm that guarantees a constrained estimate.

``` r
Integrating piecewise results 
Estimated rg is within (-1,1). Continuing computing standard error. 
```

The last section gives the genetic correlation, its standard error and P-value based on Wald test.

``` r
Genetic Correlation:  -0.2354 (0.0444)
P:  1.146e-07
```

For Help
--------

If you have specific questions, you may email the maintainer of `HDL` via <zheng.ning@ki.se>.
