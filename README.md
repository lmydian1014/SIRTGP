# SIRTGP

An R package for performing the Bayesian time-varying classification with signal interactions via relaxed thresholded Gaussian Process
### Install and load the package
```
list.of.packages = c("truncnorm", "MASS", "invgamma", "BayesGPfit", "lattice", "sn", "coda", "mcmcplots", "mcmc","Rcpp", "plgp", "matrixStats", "igraph", "easypackages")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
library(easypackages)
libraries(list.of.packages)

devtools::install_github("lmydian1014/SIRTGP")
library(SIRTGP)
