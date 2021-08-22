HDL has been updated to v1.4.0 (2021-04-15)
------------

* Add LD reference building pipline. Please see [Build a reference panel](https://github.com/zhenin/HDL/wiki/Build-a-reference-panel) for more information.
* Adjust code of loading LD reference for the compatibility of customized LD reference.

HDL has been updated to v1.3.10 (2020-04-05)
------------
* A new optional argument `fill.missing.N` is added in both the R version and the command line version. If NULL (default), the SNPs with missing N are removed. One of "median", "min" or "max" can be given so that the missing N will be filled accordingly. For example, "median" means the missing N are filled with the median N of the SNPs with available N.
* A new optional argument `intercept.output` is added in the R version for those who are interested in. By setting `intercept.output=TRUE`, the intercept terms estimates and se will be reported in estimates.df.
* A bug related to multiallelic or duplicated SNPs is fixed.
* Codes have been added to handle when column `b` is OR instead of log(OR).

HDL has been updated to v1.3.9 (2020-11-24)
------------
*   SNPs without N information are now excluded during data loading. 
*   A bug related to lower case A1 and A2 is fixed. The bug might lead to h2 mistakenly estimated to be 0.


HDL has been updated to v1.3.8 (2020-07-12)
------------
* A new optional argument `jackknife.df` is added. By setting `jackknife.df=TRUE`, it allows to output the block-jackknife estimates of heritabilities, genetic covariance and genetic correlation. 
*   Please see [here](https://github.com/zhenin/HDL/wiki/Syntax-and-results-of-HDL#Estimating-genetic-correlation-using-HDL) for more details about the new argument `jackknife.df`.
*   Full change log can be found [here](https://github.com/zhenin/HDL/blob/master/HDL/ChangeLog).


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

HDL in a nut shell
-----------------------
A short presentation about the main ideas and results of HDL given at EMGM 2020 is available [here](https://www.youtube.com/watch?v=Q2VR1iL4l9o&list=PLnJh2XY-rMTnOdlGMEUgJfoS_Uu9Qcun0&index=6&t=0s) from 37 to 49 minutes.


Citation
--------
If you use the HDL software, please cite

[Ning, Z., Pawitan, Y. & Shen, X. High-definition likelihood inference of genetic correlations across human complex traits. _Nat Genet_ (2020)](https://www.nature.com/articles/s41588-020-0653-y).


For Help
--------

For direct R documentation of `HDL.rg` function, you can use a question mark in R:

``` r
?HDL.rg
```

Some bugs might have been reported and solved in the latest version of `HDL`. Therefore, please make sure your `HDL` has been updated to the latest version (see [here](https://github.com/zhenin/HDL/wiki/Installation-and-update#Updating-HDL) for how to update `HDL`).

If you have questions, you may find the [FAQ](https://github.com/zhenin/HDL/wiki/FAQ) page is helpful. If you want further discussion or still have questions, please feel free to email the maintainer of `HDL` via <xia.shen@ki.se>.

Acknowledgement
--------

Thank all of you who have supported this project or reported bugs! Special thanks to Dr. Paul RHJ Timmers (The University of Edinburgh) for his active bug reporting of early HDL versions.
