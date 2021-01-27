For reference, this directory contains the expected output from TSexample.  It reproduces the run used in the paper (with the packages that were current at the time). 

"example-output.Rdata" contains: 
  all.meta.out = all.meta with the TS prediction column obtained in the run;
  session = the session info for the run.

"Fig1-output.pdf" is the version of Fig 1 generated in the run (with default R colors).

----------------------------------------------------------------------

Note that this also works with updated versions of R, as illustrated with the session info below:

> source("TSexample.R") # rerun on a new installation
Loading required package: Matrix
Loaded glmnet 4.1
Warning message:
In RNGkind("Mersenne-Twister", "Inversion", "Rounding") :
  non-uniform 'Rounding' sampler used

> load("example-output/example-output.Rdata") # load original run

> identical(all.meta,all.meta.out) # identical results?
[1] TRUE

> session # ORIGINAL session info
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

Random number generation:
 RNG:     
 Normal:  
 Sample:  
 
locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] glmnet_2.0-16  foreach_1.4.4  Matrix_1.2-18  limma_3.36.5   colorout_1.0-2

loaded via a namespace (and not attached):
[1] compiler_3.5.1   tools_3.5.1      codetools_0.2-15 grid_3.5.1      
[5] iterators_1.0.10 lattice_0.20-35 

> sessionInfo() # NEW session info
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] glmnet_4.1     Matrix_1.3-0   limma_3.44.3   colorout_1.0-2

loaded via a namespace (and not attached):
 [1] compiler_4.0.3   tools_4.0.3      survival_3.2-7   splines_4.0.3   
 [5] codetools_0.2-18 grid_4.0.3       iterators_1.0.13 foreach_1.5.1   
 [9] shape_1.4.5      lattice_0.20-41 

