
######## generate data from BCI true model
gen_data_sim = function(K = 6, rt = 30, num_charater = 19, train_seq = 15, test_seq = 5, sd = 5){
    beta_tar = matrix(0, nrow = K, ncol = rt)
    beta_ntar = matrix(0, nrow = K, ncol = rt)

    Sigma_tar =  matrix(c(1, 0.7, 0.1, 0.1, 0.7, 1,0.1,0.1,0.1,0.1,1,0.7,0.1,0.1,0.7,1),K,K)
    Sigma_ntar = matrix(c(1, 0.1, 0.1, 0.5, 0.1, 1,0.1,0.1,0.1,0.1,1,0.1,0.5,0.1,0.1,1),K,K)
     
    Y_train = rep(c(1,rep(0, 5), 1, rep(0, 5)), train_seq*num_charater)
    Y_test = rep(c(1,rep(0, 5), 1, rep(0, 5)), test_seq*num_charater)
    charater_train = sort(rep(c(0:18), 12*train_seq))
    charater_test = sort(rep(c(0:18), 12*test_seq))

    code_train = rep(c(1:12), train_seq * num_charater)
    code_test = rep(c(1:12), test_seq * num_charater)

    beta_tar[1,] = c(0.000000e+00, 8.229730e-01, 1.623497e+00, 2.379737e+00, 3.071064e+00, 3.678620e+00, 4.185832e+00, 4.578867e+00, 4.847001e+00, 
                    4.982922e+00, 4.982922e+00, 4.847001e+00, 4.578867e+00, 4.185832e+00, 3.678620e+00, 3.071064e+00, 2.379737e+00,  1.623497e+00,  
                    8.229730e-01,  6.123234e-16, -5.877853e-01, -9.510565e-01, -9.510565e-01, -5.877853e-01, -2.449294e-16,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    
    beta_ntar[1,] = c(0.000000e+00,  1.645946e-01,  3.246995e-01,  4.759474e-01,  6.142127e-01,
                 7.357239e-01,  8.371665e-01,  9.157733e-01,  9.694003e-01,  9.965845e-01,
                 9.965845e-01,  9.694003e-01,  9.157733e-01, 8.371665e-01,  7.357239e-01,
                 6.142127e-01,  4.759474e-01,  3.246995e-01,  1.645946e-01,  1.224647e-16,
                -1.175571e-01, -1.902113e-01, -1.902113e-01, -1.175571e-01, -4.898587e-17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    beta_tar[2,] = c(0, 0, 0, 0, 0, -2.449294e-16, -5.877853e-01, -9.510565e-01, -9.510565e-01,
                 -5.877853e-01,  6.123234e-16,  8.229730e-01,  1.623497e+00,
                 2.379737e+00,  3.071064e+00,  3.678620e+00,  4.185832e+00,
                 4.578867e+00,  4.847001e+00,  4.982922e+00,  4.982922e+00,
                 4.847001e+00,  4.578867e+00,  4.185832e+00,  3.678620e+00,
                 3.071064e+00,  2.379737e+00,  1.623497e+00,  8.229730e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    beta_ntar[2,] = c(0, 0, 0, 0, 0, -4.898587e-17, -1.175571e-01, -1.902113e-01, -1.902113e-01,
                 -1.175571e-01,  1.224647e-16,  1.645946e-01,  3.246995e-01,
                 4.759474e-01,  6.142127e-01,  7.357239e-01,  8.371665e-01,
                 9.157733e-01,  9.694003e-01,  9.965845e-01,  9.965845e-01,
                 9.694003e-01,  9.157733e-01,  8.371665e-01,  7.357239e-01,
                 6.142127e-01,  4.759474e-01,  3.246995e-01,  1.645946e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    beta_tar[3,] = c(0, 0, 0, 0, 0,0, 0, 0, 0.000000e+00,  8.229730e-01,  1.623497e+00,  2.379737e+00,  3.071064e+00,
                 3.678620e+00,  4.185832e+00,  4.578867e+00,  4.847001e+00,  4.982922e+00,
                 4.982922e+00,  4.847001e+00,  4.578867e+00,  4.185832e+00,  3.678620e+00,
                 3.071064e+00,  2.379737e+00,  1.623497e+00,  8.229730e-01,  6.123234e-16,
                 -5.877853e-01, -9.510565e-01, -9.510565e-01, -5.877853e-01, -2.449294e-16, 0, 0, 0, 0, 0, 0, 0)

    beta_ntar[3,] = c(0, 0, 0, 0, 0,0, 0, 0, 0.000000e+00,  1.645946e-01,  3.246995e-01,  4.759474e-01,  6.142127e-01,
                 7.357239e-01,  8.371665e-01,  9.157733e-01,  9.694003e-01,  9.965845e-01,
                 9.965845e-01,  9.694003e-01,  9.157733e-01, 8.371665e-01,  7.357239e-01,
                 6.142127e-01,  4.759474e-01,  3.246995e-01,  1.645946e-01,  1.224647e-16,
                -1.175571e-01, -1.902113e-01, -1.902113e-01, -1.175571e-01, -4.898587e-17, 0, 0, 0, 0, 0, 0, 0)


    beta_tar[4,] = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.449294e-16, -5.877853e-01, -9.510565e-01, -9.510565e-01,
                 -5.877853e-01,  6.123234e-16,  8.229730e-01,  1.623497e+00,
                 2.379737e+00,  3.071064e+00,  3.678620e+00,  4.185832e+00,
                 4.578867e+00,  4.847001e+00,  4.982922e+00,  4.982922e+00,
                 4.847001e+00,  4.578867e+00,  4.185832e+00,  3.678620e+00,
                 3.071064e+00,  2.379737e+00,  1.623497e+00,  8.229730e-01, 0,  0,  0)
    beta_ntar[4,] = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0, -4.898587e-17, -1.175571e-01, -1.902113e-01, -1.902113e-01,
                 -1.175571e-01,  1.224647e-16,  1.645946e-01,  3.246995e-01,
                 4.759474e-01,  6.142127e-01,  7.357239e-01,  8.371665e-01,
                 9.157733e-01,  9.694003e-01,  9.965845e-01,  9.965845e-01,
                 9.694003e-01,  9.157733e-01,  8.371665e-01,  7.357239e-01,
                 6.142127e-01,  4.759474e-01,  3.246995e-01,  1.645946e-01, 0,  0,  0)
    # beta_tar[5,] = beta_tar[2,]/4
    # beta_ntar[5,] = beta_tar[2,]/4
    # #beta_tar[6,] = beta_tar[4,]/4
    # #beta_ntar[6,] = -beta_ntar[4,]
    # beta_tar[6,] = c(0, 0, 0, 0, 0, 0, 0, 0, 0, -2.449294e-16, -5.877853e-01, -9.510565e-01, -9.510565e-01,
    #              -5.877853e-01,  6.123234e-16,  8.229730e-01,  1.623497e+00,
    #              2.379737e+00,  3.071064e+00,  3.678620e+00,  4.185832e+00,
    #              4.578867e+00,  4.847001e+00,  4.982922e+00,  4.982922e+00,
    #              4.847001e+00,  4.578867e+00,  4.185832e+00,  3.678620e+00,
    #              3.071064e+00,  2.379737e+00,  1.623497e+00,  8.229730e-01, 0,  0,  0, 0, 0, 0, 0, )
    # beta_ntar[6,] = -c(0, 0, 0, 0, 0, 0, 0,0,0, -4.898587e-17, -1.175571e-01, -1.902113e-01, -1.902113e-01,
    #              -1.175571e-01,  1.224647e-16,  1.645946e-01,  3.246995e-01,
    #              4.759474e-01,  6.142127e-01,  7.357239e-01,  8.371665e-01,
    #              9.157733e-01,  9.694003e-01,  9.965845e-01,  9.965845e-01,
    #              9.694003e-01,  9.157733e-01,  8.371665e-01,  7.357239e-01,
    #              6.142127e-01,  4.759474e-01,  3.246995e-01,  1.645946e-01, 0,  0,  0,0, 0, 0, 0,)
    ##  beta_tar is 30*4
    beta_tar = beta_tar/4
    X_train = matrix(0, nrow = length(Y_train), ncol = rt*K) ##3420 (12*15*19)  120 (30*4)
    X_test = matrix(0, nrow = length(Y_test), ncol = rt*K) ##3420 (12*15*19)  120 (30*4)

    for(i in c(1:length(Y_train))){
        eps_tar = mvrnorm(n = rt, rep(0, K), Sigma_tar)
        eps_ntar = mvrnorm(n = rt, rep(0, K), Sigma_ntar)
        eps = rnorm(rt*K, 0, sd)
        X_train[i,] = c(beta_tar) * Y_train[i] + c(beta_ntar) * (1-Y_train[i]) + c(eps_tar) * Y_train[i] + c(eps_ntar) * (1-Y_train[i]) + eps
    }
    for(i in c(1:length(Y_test))){
        eps_tar = mvrnorm(n = rt, rep(0, K), Sigma_tar)
        eps_ntar = mvrnorm(n = rt, rep(0, K), Sigma_ntar)
        eps = rnorm(rt*K, 0, sd)
        X_test[i,] = c(beta_tar) * Y_test[i] + c(beta_ntar) * (1-Y_test[i]) + c(eps_tar) * Y_test[i] + c(eps_ntar) * (1-Y_test[i]) + eps
    }
    X_train = t(X_train)
    X_test = t(X_test)

    mean_train = rowMeans(X_train)
    std_train = rowSds(X_train)
    X_train = t(standard(t(X_train)))  ### p* n
    X_test = (X_test - mean_train)/std_train  

    n = ncol(X_train)
    n_test = ncol(X_test)

    V0 = K*(K-1)/2
    X0_train = matrix(0, nrow = V0, ncol = n)
    X0_test = matrix(0, nrow = V0, ncol = n_test)
    for(i in c(1:n)){
       for(u in c(1:(K-1))){
          for(v in c((u+1):K)){
             X0_train[(2*K-2-u)*(u-1)/2 + v-1,i] = cor(X_train[(rt*(u-1)+1):(rt*u),i], X_train[(rt*(v-1)+1):(rt*v),i])
          }
       }
    }

    for(i in c(1:n_test)){
       for(u in c(1:(K-1))){
          for(v in c((u+1):K)){
             X0_test[(2*K-2-u)*(u-1)/2 + v-1,i] = cor(X_test[(rt*(u-1)+1):(rt*u),i], X_test[(rt*(v-1)+1):(rt*v),i])
          }
       }
    }

    X0_train_F = FisherZ(X0_train)
    X0_test_F = FisherZ(X0_test)
    mean_train0 = rowMeans(X0_train_F)
    std_train0 = rowSds(X0_train_F)

    X0_train_final = t(standard(t(X0_train_F))) ###  120 3420 V0*n

    X0_test_final = (X0_test_F - mean_train0)/std_train0  ####  
    X0_train = X0_train_final  ## 6 3420
    X0_test = X0_test_final ##6 1140
    return(list('X_train' = X_train, 'X_test' = X_test, 'X0_train' = X0_train, 'X0_test' = X0_test, 'Y_train' = Y_train, 'Y_test' = Y_test, 
        'charater_train' = charater_train, 'charater_test' = charater_test, 'code_train' = code_train, 'code_test' = code_test))
}




gen_data = function(n_train = 100, n_test = 100, d = 1L, K = 16, rt = 26, grids_lim = c(-1,1), random = FALSE, poly_degree = 8L,
              a = 0.1, b = 1, center = NULL, rate = NULL, max_range = 6, thres1 = 0, thres2 = 0, tausq = 0.01, sigsq = 0.0001){
    # 1 for train 2 for test
    set.seed(1)
    ############### spatial-wise ###############
    grids = GP.generate.grids(d = d, num_grids = rt, grids_lim = grids_lim)
    V = nrow(grids)
    lambda = GP.eigen.value(poly_degree=poly_degree, a=a, b=b, d=d)
    Xmat = GP.eigen.funcs.fast(grids = grids, poly_degree=poly_degree, a=a, b=b)
    L = length(lambda)
    V0 = K*(K-1)/2
    ############### time-wise ###############
    e = matrix(nrow = L, ncol = K)
    E = matrix(nrow = rt, ncol = K)
    E_hat = matrix(nrow = rt, ncol = K)
    beta = matrix(nrow = rt, ncol = K)
    
    eta = rnorm(V0, 0, 1)
    eta_hat = rnorm(V0, mean = eta, sd = sqrt(sigsq))
    for(iter in c(1:K)){
        e[,iter] = rnorm(L, 0, sqrt(lambda))
    }
    E = Xmat %*% e ### (T*L)*(L*K) = T*K
    
    E_hat = matrix(rnorm(length(E), mean = E, sd = sqrt(sigsq)), ncol = K)
    
    beta = E/(K*rt) * G_thres(abs(E_hat), thres1)
    beta0 = eta/V0 * G_thres(abs(eta_hat), thres2)  

    X_train = matrix(nrow = rt*K, ncol = n_train)
    X_test = matrix(nrow = rt*K, ncol = n_test)
    for(iter in c(1:K)){
        print(iter)
        X_train[((iter-1)*rt+1) : (iter*rt),] = Xmat %*% matrix(rnorm(n_train*L,mean=0,sd=sqrt(lambda)), nrow = L) ## (T*L) * (L*n)
        X_test[((iter-1)*rt+1) : (iter*rt),] = Xmat %*% matrix(rnorm(n_test*L,mean=0,sd=sqrt(lambda)), nrow = L) ## (T*L) * (L*n_test)
    }

    X0_train = matrix(NA, nrow = K*(K-1)/2, ncol = n_train)
    X0_test = matrix(NA, nrow = K*(K-1)/2, ncol = n_test)
   
    for(i in c(1:n_train)){
        for(u in c(1:(K-1))){
            for(v in c((u+1): K)){
                id = (2*K-2-u)*(u-1)/2 + v-1
                X0_train[id,i] = cor(X_train[((u-1)*rt+1):(u*rt),i], X_train[((v-1)*rt+1):(v*rt),i])
            }
        }
    }
    for(i in c(1:n_test)){
        for(u in c(1:(K-1))){
            for(v in c((u+1): K)){
                id = (2*K-2-u)*(u-1)/2 + v-1
                X0_test[id,i] = cor(X_test[((u-1)*rt+1):(u*rt),i], X_test[((v-1)*rt+1):(v*rt),i])
            }
        }
    }
    Y_train = rep(0, n_train)
    Y_test = rep(0, n_test)
    
    for(iter in c(1:K)){
        Y_train = Y_train + colSums(c(beta[,iter]) * X_train[((iter-1)*rt+1):(iter*rt), ])
        Y_test = Y_test + colSums(c(beta[,iter]) * X_test[((iter-1)*rt+1):(iter*rt), ])
    }
    Y_train = Y_train + colSums(beta0 * X0_train)
    Y_test = Y_test + colSums(beta0 * X0_test)
    #Y = Y + (beta0 %*% X0)
    eps = rnorm(n_train, mean = 0, sd = sqrt(tausq))
    Y_train = Y_train + eps

    E_hat_v = c()
    for(i in c(1:ncol(E_hat))){
        E_hat_v = c(E_hat_v, E_hat[,i])
    }
    p = pnorm(Y_train-eps)
    # print(range(Y))
    # print(range(eps))
    # print(range(p))

    #R_sq = sum((p - mean(p))^2)/sum(p*(1-p))
    #SNR = R_sq/(1-R_sq)
    
    SNR = sum((Y_train-eps - mean(Y_train))^2)/sum(eps^2)
    R_sq = SNR/(1+SNR)
    Y_train_cat = Y_train
    Y_test_cat = Y_test
    for(i in c(1:n_train)){
        if(Y_train[i]>0){
            Y_train_cat[i] = 1
        }
        else{
            Y_train_cat[i] = 0
        }
    }
    for(i in c(1:n_test)){
        if(Y_test[i]>0){
            Y_test_cat[i] = 1
        }
        else{
            Y_test_cat[i] = 0
        }
    }

    return(list('grids' = grids, 'X_train' = X_train, 'X0_train' = X0_train, 'X_test' = X_test, 'X0_test' = X0_test, 'Y_train' = Y_train_cat, 'Y_test' = Y_test_cat, 
        'lambda' = lambda, 'Xmat' = Xmat, 'E' = E, 'E_hat' = E_hat_v, 'eta' = eta, 'eta_hat' = eta_hat, 
      'e' = e, 'beta' = beta, 'beta0' = beta0, 'SNR' = SNR, 'R_sq' = R_sq, 'Z' = Y_train-eps, 'eps' = eps,'p' = p))
}




G_thres = function(x, thres){
	return(ifelse(x > thres, 1, 0))
}
 





