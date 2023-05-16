

analysis_chain_bci_no_M = function(T, burn_in, K, rt, n, chain, X, X0, Xmat){
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


analysis_chain_bci_v3 = function(T, burn_in, K, n, chain, X, X0, Xmat){
    #e_record = chain$e
    rt = 128
    e_record = list()
    e_record[[1]] = chain$e[,1:K]
    E_hat = chain$E_hat
    eta_hat = chain$eta_hat
    thres1 = chain$thres1
    thres2 = chain$thres2
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

        eta_m_temp = matrix(0, nrow = K, ncol = K)
        #eta_m_temp[lower.tri(eta_m_temp, diag=FALSE)] = eta_temp
        eta_m_temp[lower.tri(eta_m_temp, diag=FALSE)] = eta_temp  * G_thres(abs(eta_hat[,t]), thres2[t])
        eta_m_temp = t(eta_m_temp)
        eta_m_temp = eta_m_temp + t(eta_m_temp)
        M = rowSums(eta_m_temp)
        
        R = list()
        E = list()
        beta = list()
        for(iter in c(1:K)){
            R[[iter]] = Xmat %*% e_record[[t]][,iter]
            E[[iter]] = M[iter] + R[[iter]]
            beta[[iter]] = E[[iter]]/(K*rt) * G_thres(abs(E_hat[((iter-1)*rt+1):(iter*rt),t]), thres1[t])
        }
        beta0 = eta_temp/V0 * G_thres(abs(eta_hat[,t]), thres2[t])  ## with length 6

        
        for(iter in c(1:K)){
            Y[,t] = Y[,t] + colSums(c(beta[[iter]]) * X[((iter-1)*rt+1):(iter*rt),])
        }
        Y[,t] = Y[,t] + colSums(beta0 * X0)
        beta0_record[,t] = beta0
        beta_record[,t] = unlist(beta)

    }
    return(list(Y, beta0_record, beta_record))
}

analysis_chain_bci_nothres = function(T, burn_in, K, n, chain, X, X0, Xmat){
    #e_record = chain$e
    rt = 26
    e_record = list()
    e_record[[1]] = chain$e[,1:K]
    E_hat = chain$E_hat
    eta_hat = chain$eta_hat
    

    for(t in c(2:T)){
      #print(t)
      e_record[[t]] = chain$e[,((t-1)*K+1) : (t*K)]
    }
    
    Y = matrix(0, nrow = n, ncol = T)
    for(t in c((burn_in + 1):T)){
        eta_temp = chain$eta[,t]

        eta_m_temp = matrix(0, nrow = K, ncol = K)
        eta_m_temp[lower.tri(eta_m_temp, diag=FALSE)] = eta_temp
        eta_m_temp = t(eta_m_temp)
        eta_m_temp = eta_m_temp + t(eta_m_temp)
        M = rowSums(eta_m_temp)
        
        R = list()
        E = list()
        beta = list()
        for(iter in c(1:K)){
            R[[iter]] = Xmat %*% e_record[[t]][,iter]
            E[[iter]] = M[iter] + R[[iter]]
            beta[[iter]] = E[[iter]] * G_thres(abs(E_hat[((iter-1)*rt+1):(iter*rt),t]), thres)
        }
        beta0 = eta_temp * G_thres(abs(eta_hat[,t]), thres)  ## with length 6

        
        for(iter in c(1:K)){
            Y[,t] = Y[,t] + colSums(c(beta[[iter]]) * X[((iter-1)*rt+1):(iter*rt),])
        }
        Y[,t] = Y[,t] + colSums(beta0 * X0)

    }
    return(list(Y, beta0, beta))
}


analysis_chain_bci = function(T, burn_in, K, n, chain, X, X0, Xmat){
    #e_record = chain$e
    rt = 26
    e_record = list()
    e_record[[1]] = chain$e[,1:K]
    E_hat = chain$E_hat
    eta_hat = chain$eta_hat
    thres1 = chain$thres1
    thres2 = chain$thres2

    for(t in c(2:T)){
      #print(t)
      e_record[[t]] = chain$e[,((t-1)*K+1) : (t*K)]
    }
    
    Y = matrix(0, nrow = n, ncol = T)
    for(t in c((burn_in + 1):T)){
        eta_temp = chain$eta[,t]

        eta_m_temp = matrix(0, nrow = K, ncol = K)
        eta_m_temp[lower.tri(eta_m_temp, diag=FALSE)] = eta_temp
        eta_m_temp = t(eta_m_temp)
        eta_m_temp = eta_m_temp + t(eta_m_temp)
        M = rowSums(eta_m_temp)
        
        R = list()
        E = list()
        beta = list()
        for(iter in c(1:K)){
            R[[iter]] = Xmat %*% e_record[[t]][,iter]
            E[[iter]] = M[iter] + R[[iter]]
            beta[[iter]] = E[[iter]] * G_thres(abs(E_hat[((iter-1)*rt+1):(iter*rt),t]), thres1[t])
        }
        beta0 = eta_temp * G_thres(abs(eta_hat[,t]), thres2[t])  ## with length 6

        
        for(iter in c(1:K)){

            Y[,t] = Y[,t] + colSums(c(beta[[iter]]) * X[((iter-1)*rt+1):(iter*rt),])
        }
        Y[,t] = Y[,t] + colSums(beta0 * X0)

    }
    return(list(Y, beta0, beta))
}




analysis_chain = function(T, K, n, chain, X, X0, Xmat, Xmat0, region){
    g_record = chain$g
    #e_record = chain$e
    e_record = list()
    e_record[[1]] = chain$e[,1:K]
    for(t in c(2:T)){
      #print(t)
      e_record[[t]] = chain$e[,((t-1)*K+1) : (t*K)]
    }
    
    rc = c()
    for(iter in c(1:K)){
        rc = c(rc, nrow(Xmat[region==iter, ]))
    }
    
    E_hat_record = list()
    for(t in c(1:T)){
        E_hat_record[[t]] = list()
        for(iter in c(1:K)){
            if(iter == 1){
                ind_E_hat = seq(1, rc[1], by = 1);   
            }
            else{
                ind_E_hat = seq(sum(rc[1: (iter-1)]) + 1, sum(rc[1:iter]), by = 1);      
            }  
            E_hat_record[[t]][[iter]] = c(chain$E_hat[ind_E_hat, t])
        }
    }


    #E_hat_record = chain$E_hat
    eta_hat_record = chain$eta_hat
    L = nrow(g_record)
    Y = matrix(0, nrow = n, ncol = T)
    for(t in c(1:T)){
        eta = c()
        M = rep(0, K)
        beta = list()
        for(u in c(1:(K-1))){
            for(v in c((u+1):K)){
                eta = append(eta, (Xmat0[u,] * Xmat0[v,]) %*% g_record[,t])
            }
        }


        # for(u in c(1:K)){
        #     for(v in c(1:K)){

        #         if(v < u){
        #             index = (2*K-2-v)*(v-1)/2 + u-1
        #             M[u] = M[u] + eta[index]
        #         }
        #         if(u < v){ 
        #             index = (2*K-2-u)*(u-1)/2 + v-1
        #             #print(u)
        #             #print(v)
        #             M[u] = M[u] + eta[index]
        #         }
        #     }
        # }
        R = list()
        E = list()
        for(iter in c(1:K)){
            M[iter] = sum(colSums(Xmat0[-iter, ]) * Xmat0[iter, ] * g_record[,t])
            R[[iter]] = Xmat[region==iter, ] %*% e_record[[t]][,iter]
            E[[iter]] = M[iter] + R[[iter]]
            beta[[iter]] = E[[iter]] * G_thres(abs(E_hat_record[[t]][[iter]]), thres)
        }
        beta0 = eta * G_thres(abs(eta_hat_record[,t]), thres)  ## with length 6

        
        for(k in c(1:K)){
            Y[,t] = Y[,t] + colSums(c(beta[[k]]) * X[region==k,])
        }
        Y[,t] = Y[,t] + colSums(beta0 * X0)

    }
    return(Y)
}


# n_test = 100
# X_test = matrix(NA, nrow = nrow(grids), ncol = n_test)
# X0_test = matrix(NA, nrow = K*(K-1)/2, ncol = n_test)
# X_test = Xmat %*% matrix(rnorm(n_test*L,mean=0,sd=sqrt(lambda)), nrow = L)


# for(i in c(1:n_test)){
#     for(u in c(1:(K-1))){
#         for(v in c((u+1): K)){
#             id = (2*K-2-u)*(u-1)/2 + v-1
#             X0_test[id,i] = cor(X_test[region == u,i], X_test[region == v,i])
#         }
#     }
# }

# Y_test = rep(0, n)
#  for(iter in c(1:K)){
#     Y_test = Y_test + colSums(c(dat$beta[[iter]]) * X_test[region==iter,])
# }
# Y_test = Y_test + colSums(dat$beta0 * X0_test)
# #Y = Y + (beta0 %*% X0)
# eps = rnorm(n, mean = 0, sd = sqrt(tausq))
# Y = Y + eps






analysis_chain_v3 = function(T, K, n, chain, X, X0, Xmat, Xmat0, region, E_hat, eta_hat){
    #e_record = chain$e
    e_record = list()
    E_hat_record = list()
    e_record[[1]] = chain$e[,1:K]
    E_hat_record[[1]] = chain$E_hat[,1]

    for(t in c(2:T)){
      #print(t)
      e_record[[t]] = chain$e[,((t-1)*K+1) : (t*K)]
      E_hat_record[[iter]] = E_hat[((iter - 1)*rt + 1): (iter*rt)]
    }
    
    E_hat_record = list()
    for(iter in c(1:K)){
        E_hat_record[[iter]] = E_hat[((iter - 1)*rt + 1): (iter*rt)]
    }

    L = nrow(e_record)
    Y = matrix(0, nrow = n, ncol = T)
    for(t in c(1:T)){
        eta = c()
        M = rep(0, K)
        beta = list()
        for(u in c(1:(K-1))){
            for(v in c((u+1):K)){
                eta = append(eta, (Xmat0[u,] * Xmat0[v,]) %*% g_record[,t])
            }
        }
        R = list()
        E = list()
        for(iter in c(1:K)){
            M[iter] = sum(colSums(Xmat0[-iter, ]) * Xmat0[iter, ] * g_record[,t])
            R[[iter]] = Xmat[region==iter, ] %*% e_record[[t]][,iter]
            E[[iter]] = M[iter] + R[[iter]]
            beta[[iter]] = E[[iter]] * G_thres(abs(E_hat[[iter]]), thres)
        }
        beta0 = eta * G_thres(abs(eta_hat), thres)  ## with length 6

        
        for(k in c(1:K)){
            Y[,t] = Y[,t] + colSums(c(beta[[k]]) * X[region==k,])
        }
        Y[,t] = Y[,t] + colSums(beta0 * X0)

    }
    return(Y)
}




