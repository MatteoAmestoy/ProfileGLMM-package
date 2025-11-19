## R CMD check results

0 errors | 0 warnings | 1 note

*1. Purpose

The package implements profile regression estimation for Generalized Linear Mixed Models (GLMMs).


*2. C++ Standard and Requirements

The package uses C++ for performance, linking against Rcpp, RcppArmadillo, and RcppDist.

*3. Non-Standard Dependencies

The package imports the following specialized MCMC sampler packages. 
These are required for internal utility functions that are part of the overall fitting procedure:

LaplacesDemon: Used for the rinvwishart sampler.

MCMCpack: Used for the rdirichlet sampler.
