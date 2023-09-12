

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




