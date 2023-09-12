\name{SIRTGP_fit}
\alias{SIRTGP_fit}
\title{
Bayesian fitting of the time-varying regression with signal interactions via the Relaxed Thresholded Gaussian Process
}

\usage{
SIRTGP_fit(T, K, L, n, Y, X, X0, Xmat, eta, eta_m, e, E_hat, vec eta_hat, double beta0,
	double thres1, double thres2, vec lambda, double tausq, double sigsq, double sigsq_eta, int rt, double alpha1, double alpha2, vec prob1, vec prob2, vec val1, vec val2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{A integer number to specify the number of iteractions in MCMC sampling.}
  \item{K}{A integer number to specify the number of channels.}
  \item{L}{A integer number to specify the number of basis functions.}
  \item{n}{A integer number to specify the number of flashes.}
  \item{Y}{A binary vector with length n representing the target/non-target flash.}
  \item{X}{A matrix representing the EEG signal with dimension K*rt by n, where rt is the number of time points collected on each channel.}
  \item{X0}{A matrix representing the EEG signal interaction with dimension K*(K-1)/2 by n.}
  \item{Xmat}{A matrix represents the basis functions evaluated at the grid points, where rows are observations and columns are the basis functions.}
  \item{eta}{A vector with length K*(K-1)/2 specifies the initial value of eta.}
  \item{eta_m}{A matrix with dimension K by K specifies the initial value of the matrix version of eta.}
  \item{e}{A matrix with dimension L by K pecifies the initial value of the matrix version of e.}
  \item{E_hat}{A vector with length K*rt specifies the initial value of E_hat.}
  \item{eta_hat}{A vector with length K*(K-1)/2 specifies the initial value of eta_hat.}
  \item{beta0}{A scalar represents the initial value of the intercept.}
  \item{thres1}{A scalar represents the initial value of the 1st thresholding value in the main effect part.}
  \item{thres2}{A scalar represents the initial value of the 2nd thresholding value in the interaction effect part.}
  \item{lambda}{A vector with length L represents the eigen values.}
  \item{tausq}{A scalar represents the variance of the noise \epsilon.}
  \item{sigsq}{A scalar represents the variance in the probit model.}
  \item{sigsq_eta}{A scalar represents the intial value of sigsq_eta.}
  \item{rt}{A integer represents the number of time points collected on each channel.}
  \item{prob1}{A vector represents the initial probability distribution for the 1st thresholding values.}
  \item{prob2}{A vector represents the initial probability distribution for the 2nd thresholding values.}
  \item{val1}{A vector represents the initial value of possible choices of the 1st thresholding values.}
  \item{val2}{A vector represents the initial value of possible choices of the 2nd thresholding values.}
}

\value{
  \item{e}{A matrix of dimension L by K*T represents the posterior samples of e for each iteration.}
  \item{eta}{A matrix of dimension K*(K-1)/2 by T represents the posterior samples of eta for each iteration.}
  \item{E_hat}{A matrix of dimension K*rt by T represents the posterior samples of E_hat for each iteration.}
  \item{eta_hat}{A matrix of dimension K*(K-1)/2 by T represents the posterior samples of E_hat for each iteration.}
  \item{thres1}{A vector with length T represents the posterior samples of thres1 for each iteration.}
  \item{thres2}{A vector with length T represents the posterior samples of thres1 for each iteration.}
  \item{intercept}{A vector with length T represents the posterior samples of the intercept for each iteration.}
}
\author{
Moyan Li <moyanli@umich.edu>
}


\examples{
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
  
  T = 10
  burn_in = 5
  chain = SIRTGP_fit(T, K, L, n_train, dat$Y_train, dat$X_train, dat$X0_train, Xmat, eta_init, eta_m_init, e_init, E_hat_init, eta_hat_init, 0, thres1, thres2, lambda, tausq, sigsq, sigsq_eta, rt, alpha1, alpha2, prob1, prob2, val1, val2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
