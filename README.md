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
```

### generate data 
```
n_train = 12*19*10
n_test = 12*19*5
K = 16 #number of channels
rt = 26 #number of time points collected on each channel
tausq = 0.001; sigsq = 0.0001
thres1 = 0.0; thres2 = 0.0
dat = gen_data(n = n_train, n_test = n_test, d = 1L, K = K, rt = rt, grids_lim = c(0,1), a = 0.1, b = 10, thres1 = thres1, thres2 = thres2, tausq = tausq, sigsq = sigsq)

```

### initialize parameters
```
Xmat = dat$Xmat
lambda = dat$lambda
L = length(lambda)
alpha1 = 1; alpha2 = 1
prob1 = c(1,0,0); prob2 = c(1,0,0)
val1 = c(0, 0.1, 1); val2 = c(0, 0.1, 1)
e_init = dat$e
eta_init = dat$eta
sigsq_eta = sigsq
eta_m_init = matrix(0, nrow = K, ncol = K)
eta_m_init[lower.tri(eta_m_init, diag=FALSE)] = eta_init
eta_m_init = t(eta_m_init)
eta_m_init = eta_m_init + t(eta_m_init)
E_hat_init = dat$E_hat
eta_hat_init = dat$eta_hat
```

### Run SIRTGP model
```
T = 1000
burn_in = 200
chain = SIRTGP_fit(T, K, L, n_train, Y_train, X_train, X0_train, Xmat, eta_init, eta_m_init, e_init, E_hat_init, eta_hat_init, 0, thres1, thres2, lambda, 
                      tausq, sigsq, sigsq_eta, rt, alpha1, alpha2, prob1, prob2, val1, val2)
```






### Analysis the MCMC chain 
```
res = analysis_chain_bci_no_M(T, burn_in, K, rt, n_test, chain, X_test, X0_test, Xmat)
```
### calculate classfication accuracy
```
p = (pnorm(res[[1]][,(burn_in+1):T]))
Y_chain = matrix(0, nrow = n_test, ncol = (T-burn_in))
for(i in c(1:nrow(p))){
    Y_chain[i,] = rbernoulli(T-burn_in, p[i,])
}
Y_test_pred = c()
for(i in c(1:n_test)){
    if(sum(Y_chain[i,]) > ((T-burn_in)/2)){
        Y_test_pred[i] = 1
    }
    else{
        Y_test_pred[i] = 0
    }
}
classify.accuracy(Y_test_pred, dat$Y_test)
```






















