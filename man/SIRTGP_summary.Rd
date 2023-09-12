\name{SIRTGP_summary}
\alias{SIRTGP_summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The summary of the posterior sampling results of the SIRTGP.
}

\usage{
SIRTGP_summary(T, burn_in, K, rt, n, chain, X, X0, Xmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{A integer number to specify the number of iteractions in MCMC sampling.}
  \item{burn_in}{An integer number to specify the burn-in number.}
  \item{K}{A integer number to specify the number of channels.}
  \item{rt}{A integer represents the number of time points collected on each channel.}
  \item{n}{A integer number to specify the number of flashes.}
  \item{chain}{A list of the posterior sampling results obtained from SIRTGP_fit().}
  \item{X}{A matrix representing the EEG signal with dimension K*rt by n, where rt is the number of time points collected on each channel.}
  \item{X_0}{A matrix representing the EEG signal interaction with dimension K*(K-1)/2 by n.}
  \item{Xmat}{A matrix represents the basis functions evaluated at the grid points, where rows are observations and columns are the basis functions.}
}

\value{
  \item{Y}{A matrix with dimension n by T represents the predicted value of Y using T different values of posterior samples.}
  \item{beta0_record}{A matrix with dimension K*(K-1)/2 by T represents the posterior samples of the intercept beta0.}
  \item{beta_record}{A matrix with dimension K*rt by T represents the posterior samples of the intercept beta.}
}

\author{
Moyan Li <moyanli@umich.edu>
}


\examples{
### Generate data

n_train = 12*19*10
n_test = 12*19*5
K = 16 #number of channels
rt = 26 #number of time points collected on each channel
tausq = 0.001; sigsq = 0.0001
thres1 = 0.0; thres2 = 0.0
dat = gen_data(n = n_train, n_test = n_test, d = 1L, K = K, rt = rt, grids_lim = c(0,1), a = 0.1, b = 10, thres1 = thres1, thres2 = thres2, tausq = tausq, sigsq = sigsq)
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

##### Run the SIRTGP model

T = 1000
burn_in = 200
chain = SIRTGP_fit(T, K, L, n_train, dat$Y_train, dat$X_train, dat$X0_train, Xmat, eta_init, eta_m_init, e_init, E_hat_init, eta_hat_init, 0, thres1, thres2, lambda, tausq, sigsq, sigsq_eta, rt, alpha1, alpha2, prob1, prob2, val1, val2)


##### Analyze the SIRTGP results 

res = SIRTGP_summary(T, burn_in, K, rt, n_test, chain, dat$X_test, dat$X0_test, Xmat)                      

}
