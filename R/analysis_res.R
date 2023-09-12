#'@title The summary of the posterior sampling results of the SIRTGP.

#'@param T A integer number to specify the number of iteractions in MCMC sampling. 
#'@param burn_in An integer number to specify the burn-in number. 
#'@param K A integer number to specify the number of channels. 
#'@param rt A integer represents the number of time points collected on each channel. 
#'@param n A integer number to specify the number of flashes. 
#'@param chain A list of the posterior sampling results obtained from SIRTGP_fit(). 
#'@param X A matrix representing the EEG signal with dimension K*rt by n, where rt is the number of time points collected on each channel. 
#'@param X_0 A matrix representing the EEG signal interaction with dimension K*(K-1)/2 by n. 
#'@param Xmat A matrix represents the basis functions evaluated at the grid points, where rows are observations and columns are the basis functions. 

#'
#'@return A list of variables including the model fitting results
#'\describe{
#'\item{Y}{A matrix with dimension n by T represents the predicted value of Y using T different values of posterior samples.}
#'\item{beta0_record}{A matrix with dimension K*(K-1)/2 by T represents the posterior samples of the intercept beta0.}
#'\item{beta_record}{A matrix with dimension K*rt by T represents the posterior samples of the intercept beta.}
#'}
#'
#'@author Moyan Li <moyanli@umich.edu>
#'
#'@examples

#' \examples{
#' ### Generate data
#' 
#' n_train = 12*19*10
#' n_test = 12*19*5
#' K = 16 #number of channels
#' rt = 26 #number of time points collected on each channel
#' tausq = 0.001; sigsq = 0.0001
#' thres1 = 0.0; thres2 = 0.0
#' dat = gen_data(n = n_train, n_test = n_test, d = 1L, K = K, rt = rt, grids_lim = c(0,1), a = 0.1, b = 10, thres1 = thres1, thres2 = thres2, tausq = tausq, sigsq = sigsq)
#' Xmat = dat$Xmat
#' lambda = dat$lambda
#' L = length(lambda)
#' alpha1 = 1; alpha2 = 1
#' prob1 = c(1,0,0); prob2 = c(1,0,0)
#' val1 = c(0, 0.1, 1); val2 = c(0, 0.1, 1)
#' e_init = dat$e
#' eta_init = dat$eta
#' sigsq_eta = sigsq
#' eta_m_init = matrix(0, nrow = K, ncol = K)
#' eta_m_init[lower.tri(eta_m_init, diag=FALSE)] = eta_init
#' eta_m_init = t(eta_m_init)
#' eta_m_init = eta_m_init + t(eta_m_init)
#' E_hat_init = dat$E_hat
#' eta_hat_init = dat$eta_hat
#' 
#' ##### Run the SIRTGP model
#' 
#' T = 1000
#' burn_in = 200
#' chain = SIRTGP_fit(T, K, L, n_train, dat$Y_train, dat$X_train, dat$X0_train, Xmat, eta_init, eta_m_init, e_init, E_hat_init, eta_hat_init, 0, thres1, thres2, lambda, tausq, sigsq, sigsq_eta, rt, alpha1, alpha2, prob1, prob2, val1, val2)
#' 
#' 
#' ##### Analyze the SIRTGP results 
#' 
#' res = SIRTGP_summary(T, burn_in, K, rt, n_test, chain, dat$X_test, dat$X0_test, Xmat)                      
#' 
#' }
#'
#'@export


SIRTGP_summary = function(T, burn_in, K, rt, n, chain, X, X0, Xmat){
    #e_record = chain$e
    e_record = list()
    e_record[[1]] = chain$e[,1:K]
    E_hat = chain$E_hat
    eta_hat = chain$eta_hat
    thres1 = chain$thres1
    thres2 = chain$thres2
    intecept = chain$intecept
    V0 = K*(K-1)/2

    for(t in c(2:T)){
      #print(t)
      e_record[[t]] = chain$e[,((t-1)*K+1) : (t*K)]
    }
    
    Y = matrix(0, nrow = n, ncol = T)
    beta0_record = matrix(0, nrow = V0, ncol = T)
    beta_record = matrix(0, nrow = K*rt, ncol = T)
    for(t in c((burn_in + 1):T)){
        eta_temp = chain$eta[,t]
        R = list()
        E = list()
        beta = list()
        for(iter in c(1:K)){
            R[[iter]] = Xmat %*% e_record[[t]][,iter]
            E[[iter]] =  R[[iter]]
            beta[[iter]] = E[[iter]]/(K*rt) * G_thres(abs(E_hat[((iter-1)*rt+1):(iter*rt),t]), thres1[t])
        }
        beta0 = eta_temp/V0 * G_thres(abs(eta_hat[,t]), thres2[t])  ## with length 6

        
        for(iter in c(1:K)){
            Y[,t] = Y[,t] + colSums(c(beta[[iter]]) * X[((iter-1)*rt+1):(iter*rt),])
        }
       # Y[,t] = Y[,t] + colSums(beta0 * X0)
        Y[,t] = Y[,t] + colSums(beta0 * X0) + intecept[t]
        beta0_record[,t] = beta0
        beta_record[,t] = unlist(beta)
    }
    return(list(Y, beta0_record, beta_record))
}




G_thres = function(x, thres){
	return(ifelse(x > thres, 1, 0))
}


