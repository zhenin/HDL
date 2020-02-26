What is HDL?
------------

High-Definition Likelihood (HDL) is a likelihood-based method for estimating genetic correlation using GWAS summary statistics. 
Compared to [LD Score regression (LDSC)](https://github.com/bulik/ldsc), It reduces the variance of a genetic correlation estimate by about 60%. 
Here, we provide an R-based computational tool `HDL` to implement our method. Although `HDL` is written in R, 
you can use it with the command line. So no worry if you are not an R user.

In [the wiki](https://github.com/zhenin/HDL/wiki), we provide a detailed tutorial for the application of `HDL` together with real examples. 

What data are required?
-----------------------

-   `gwas1.df` and `gwas2.df`, which are two datasets including GWAS summary statistics of genetic variants for two traits. [This page](https://github.com/zhenin/HDL/wiki/Format-of-summary-statistics) describes the format of summary statistics for `HDL`, and how to perform data wrangling.

*   The eigenvalues and eigenvectors of LD matrices. For the European-ancestry population, 
we have computed the LD matrices and their eigen-decomposition from 336,000 Genomic British individuals in UK Biobank. 
You can download these pre-computed reference files following the [instruction](https://github.com/zhenin/HDL/wiki/Reference-panels) 
in the [wiki](https://github.com/zhenin/HDL/wiki).


For Help
--------

For direct R documentation of `HDL.rg` function, you can use a question mark in R:

``` r
?HDL.rg
```

If you have specific questions, you may email the maintainer of `HDL` via <zheng.ning@ki.se>.
