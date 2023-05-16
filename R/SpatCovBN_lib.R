#simulate data
library(igraph)

simul.bivariate.proc = function(xgrid,n_obs=1,a = 0.25,b = 10,alpha = 0.8,tau_scale=0.5){
  if(is.matrix(xgrid)==FALSE)
    xgrid = as.matrix(xgrid)
  
  e.degree = GP.eigen.degree(alpha=alpha,a=a,b=b,d=ncol(xgrid))  
  n = e.degree$n
  Xmat = GP.eigen.funcs(xgrid=xgrid, n =n, a =a, b=b)
  lambda = GP.eigen.value(n=n,a=a,b=b,d=d)
  print(length(lambda))
  xi = Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
  log_tau_1_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  log_tau_2_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  tau_1_sq = exp(log_tau_1_sq)
  tau_2_sq = exp(log_tau_2_sq)
  thres = runif(1,0,1)
  sigma_pos_sq = ifelse(xi>thres,xi,0)
  sigma_neg_sq = ifelse(xi< -thres,-xi,0)
  sigma_pos = sqrt(sigma_pos_sq)
  sigma_neg = sqrt(sigma_neg_sq)
  eta_pos = matrix(NA,nrow=nrow(xgrid),ncol=n_obs)
  eta_neg = matrix(NA,nrow=nrow(xgrid),ncol=n_obs)
  Y_1 = matrix(NA, nrow=nrow(xgrid),ncol=n_obs)
  Y_2 = matrix(NA, nrow=nrow(xgrid),ncol=n_obs)
  for(i in 1:n_obs){
    eta_pos[,i] = as.numeric(sigma_pos)*Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
    eta_neg[,i] = as.numeric(sigma_neg)*Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
    Y_1[,i] = rnorm(length(tau_1_sq),sd=sqrt(tau_1_sq))
    Y_2[,i] = rnorm(length(tau_2_sq),sd=sqrt(tau_2_sq))
  }
  Y_1 = Y_1 + eta_pos + eta_neg
  Y_2 = Y_2 + eta_pos - eta_neg
  thres_xi = sigma_pos_sq - sigma_neg_sq
  abs_thres_xi = abs(thres_xi)
  var_Y_1 = tau_1_sq+abs_thres_xi
  var_Y_2 = tau_2_sq+abs_thres_xi
  
  rho = thres_xi/sqrt(var_Y_1*var_Y_2)
  return(list(tau_1_sq=tau_1_sq,tau_2_sq=tau_2_sq,xi=xi,lambda=lambda,
              rho = rho,var_Y_1 = var_Y_1,var_Y_2 = var_Y_2,xgrid=xgrid,
              Y_1 = Y_1,Y_2 = Y_2,eta_pos = eta_pos, eta_neg = eta_neg))  
}

simul.bivariate.proc.design = function(xgrid,n_obs=1,a = 0.25,b = 10,alpha = 0.8,
                                       tau_scale=0.5,pos_mag = 1, neg_mag = 1,
                                       pos_radius=0.2,neg_radius=0.1){
  if(is.matrix(xgrid)==FALSE)
    xgrid = as.matrix(xgrid)
  
  e.degree = GP.eigen.degree(alpha=alpha,a=a,b=b,d=ncol(xgrid))  

  n = e.degree$n
  print(n)
  Xmat = GP.eigen.funcs(xgrid=xgrid ,n =n ,a =a ,b=b)
  lambda = GP.eigen.value(n=n,a=a,b=b,d=d)
  print(length(lambda))
  log_tau_1_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  log_tau_2_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  tau_1_sq = exp(log_tau_1_sq)
  tau_2_sq = exp(log_tau_2_sq)
  
  #tau_1_sq = as.matrix(rep(exp(tau_scale),length=nrow(xgrid)))
  #tau_2_sq = as.matrix(rep(exp(tau_scale),length=nrow(xgrid)))
  
  #thres = runif(1,0,1)
  #xi = Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
  #sigma_pos_sq = ifelse(xi>thres,xi,0)
  #sigma_neg_sq = ifelse(xi< -thres,-xi,0)
  
  sigma_pos_sq = ifelse(abs(xgrid[,1]-0.3)+abs(xgrid[,2]-0.7)<pos_radius,pos_mag,0)
  sigma_pos_sq = ifelse(abs(xgrid[,1]-0.7)+abs(xgrid[,2]-0.7)<pos_radius,pos_mag,sigma_pos_sq)
  sigma_pos_sq = ifelse(abs(xgrid[,1]-0.3)+abs(xgrid[,2]-0.3)<pos_radius,pos_mag,sigma_pos_sq)
  sigma_neg_sq = ifelse((xgrid[,1]-0.5)^2+(xgrid[,2]-0.5)^2<neg_radius^2,neg_mag,0)
  sigma_neg_sq = ifelse(abs(xgrid[,1]-0.7)+abs(xgrid[,2]-0.3)<neg_radius,neg_mag,sigma_neg_sq)
  
  
  
  sigma_pos = sqrt(sigma_pos_sq)
  sigma_neg = sqrt(sigma_neg_sq)
  eta_pos = matrix(NA,nrow=nrow(xgrid),ncol=n_obs)
  eta_neg = matrix(NA,nrow=nrow(xgrid),ncol=n_obs)
  Y_1 = matrix(NA, nrow=nrow(xgrid),ncol=n_obs)
  Y_2 = matrix(NA, nrow=nrow(xgrid),ncol=n_obs)
  for(i in 1:n_obs){
    eta_pos[,i] = as.numeric(sigma_pos)*Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
    eta_neg[,i] = as.numeric(sigma_neg)*Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
    Y_1[,i] = rnorm(length(tau_1_sq),sd=sqrt(tau_1_sq))
    Y_2[,i] = rnorm(length(tau_2_sq),sd=sqrt(tau_2_sq))
  }
  Y_1 = Y_1 + eta_pos + eta_neg
  Y_2 = Y_2 + eta_pos - eta_neg
  thres_xi = sigma_pos_sq - sigma_neg_sq
  abs_thres_xi = abs(thres_xi)
  var_Y_1 = tau_1_sq+abs_thres_xi
  var_Y_2 = tau_2_sq+abs_thres_xi
  rho = thres_xi/sqrt(var_Y_1*var_Y_2)
  
  cor_type = ifelse(rho>0,1,0)
  cor_type = ifelse(rho<0,-1,cor_type)
  
  return(list(tau_1_sq=tau_1_sq,tau_2_sq=tau_2_sq,sigma_pos_sq = sigma_pos_sq,
              sigma_neg_sq = sigma_neg_sq,
              rho = rho,var_Y_1 = var_Y_1,var_Y_2 = var_Y_2,xgrid=xgrid,
              Y_1 = Y_1,Y_2 = Y_2,eta_pos = eta_pos, eta_neg = eta_neg,
              cor_type = cor_type))  
}


simul.bivariate.proc.design.cauchy = function(xgrid,n_obs=1,a = 0.25,b = 10,alpha = 0.8,
                                       tau_scale=0.5,pos_mag = 1, neg_mag = 1,
                                       pos_radius=0.2,neg_radius=0.1){
  if(is.matrix(xgrid)==FALSE)
    xgrid = as.matrix(xgrid)
  
  e.degree = GP.eigen.degree(alpha=alpha,a=a,b=b,d=ncol(xgrid))  
  n = e.degree$n
  Xmat = GP.eigen.funcs(xgrid=xgrid ,n =n ,a =a ,b=b)
  lambda = GP.eigen.value(n=n,a=a,b=b,d=d)
  
  log_tau_1_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  log_tau_2_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  tau_1_sq = exp(log_tau_1_sq)
  tau_2_sq = exp(log_tau_2_sq)
  
  #tau_1_sq = as.matrix(rep(exp(tau_scale),length=nrow(xgrid)))
  #tau_2_sq = as.matrix(rep(exp(tau_scale),length=nrow(xgrid)))
  
  #thres = runif(1,0,1)
  #xi = Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
  #sigma_pos_sq = ifelse(xi>thres,xi,0)
  #sigma_neg_sq = ifelse(xi< -thres,-xi,0)
  
  sigma_pos_sq = ifelse(abs(xgrid[,1]-0.3)+abs(xgrid[,2]-0.7)<pos_radius,pos_mag,0)
  sigma_pos_sq = ifelse(abs(xgrid[,1]-0.7)+abs(xgrid[,2]-0.7)<pos_radius,pos_mag,sigma_pos_sq)
  sigma_pos_sq = ifelse(abs(xgrid[,1]-0.3)+abs(xgrid[,2]-0.3)<pos_radius,pos_mag,sigma_pos_sq)
  sigma_neg_sq = ifelse((xgrid[,1]-0.5)^2+(xgrid[,2]-0.5)^2<neg_radius^2,neg_mag,0)
  sigma_neg_sq = ifelse(abs(xgrid[,1]-0.7)+abs(xgrid[,2]-0.3)<neg_radius,neg_mag,sigma_neg_sq)
  
  
  
  sigma_pos = sqrt(sigma_pos_sq)
  sigma_neg = sqrt(sigma_neg_sq)
  eta_pos = matrix(NA,nrow=nrow(xgrid),ncol=n_obs)
  eta_neg = matrix(NA,nrow=nrow(xgrid),ncol=n_obs)
  Y_1 = matrix(NA, nrow=nrow(xgrid),ncol=n_obs)
  Y_2 = matrix(NA, nrow=nrow(xgrid),ncol=n_obs)
  for(i in 1:n_obs){
    eta_pos[,i] = as.numeric(sigma_pos)*Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
    eta_neg[,i] = as.numeric(sigma_neg)*Xmat%*%rnorm(length(lambda),mean=0,sd=sqrt(lambda))
    Y_1[,i] = rnorm(length(tau_1_sq),sd=sqrt(tau_1_sq))
    Y_2[,i] = rnorm(length(tau_2_sq),sd=sqrt(tau_2_sq))
  }
  Y_1 = Y_1 + eta_pos + eta_neg + matrix(rcauchy(length(Y_1),scale = 0.0001),nrow=nrow(xgrid),ncol=n_obs)
  Y_2 = Y_2 + eta_pos - eta_neg + matrix(rcauchy(length(Y_2),scale = 0.0001),nrow=nrow(xgrid),ncol=n_obs)
  thres_xi = sigma_pos_sq - sigma_neg_sq
  abs_thres_xi = abs(thres_xi)
  var_Y_1 = tau_1_sq+abs_thres_xi
  var_Y_2 = tau_2_sq+abs_thres_xi
  rho = thres_xi/sqrt(var_Y_1*var_Y_2)
  
  cor_type = ifelse(rho>0,1,0)
  cor_type = ifelse(rho<0,-1,cor_type)
  
  return(list(tau_1_sq=tau_1_sq,tau_2_sq=tau_2_sq,sigma_pos_sq = sigma_pos_sq,
              sigma_neg_sq = sigma_neg_sq,
              rho = rho,var_Y_1 = var_Y_1,var_Y_2 = var_Y_2,xgrid=xgrid,
              Y_1 = Y_1,Y_2 = Y_2,eta_pos = eta_pos, eta_neg = eta_neg,
              cor_type = cor_type))  
}

show.simul.bivariate.proc.2D = function(sim){
  threefigs.levelplot(sim$rho,sim$var_Y_1,sim$var_Y_2,xgrid[,1],xgrid[,2],titles=c("rho","Var 1","Var 2"))
}

voxel.cor.est = function(Y_1, Y_2){
  sapply(1:nrow(Y_1),function(i) cor(Y_1[i,],Y_2[i,]))
}

voxel.cor.test = function(Y_1, Y_2){
  sapply(1:nrow(Y_1),function(i) cor.test(Y_1[i,],Y_2[i,])$p.value)
}

filter.idx = function(idx,n){
  return(ifelse(idx<=0 | idx>n,NA,idx))
}

find.neighbor.2D.eq.sp = function(num_grids,rad_voxel = 3){
  n = num_grids*num_grids
  temp = matrix(1:n,nrow=num_grids,ncol=num_grids)
  row_idx = c(row(temp))
  col_idx = c(col(temp))
  adjust_num = seq(-rad_voxel,rad_voxel,by=1)
  n_idx = expand.grid(adjust_num,adjust_num)
  n_idx = n_idx[-which(n_idx[,1]==0 & n_idx[,2]==0),]
  neighbors = NULL
  for(i in 1:nrow(n_idx)){
    neighbors = cbind(neighbors,temp[filter.idx(row_idx+n_idx[i,1]+num_grids*(col_idx-1+n_idx[i,2]),n)])
  }
  return(neighbors)
}

spatial.kernel.smooth.2D.eq.sp = function(val,xgrid,neighbors=NULL,rad_voxel=3,n_cor=0.9,cor_diff=Inf){
  if(is.null(neighbors)){
    neighbors=find.neighbor.2D.eq.sp(sqrt(nrow(xgrid)),rad_voxel)
  }
  rho = -log(n_cor)/sqrt(sum((xgrid[1,]-xgrid[2,])^2))
  neighbors = cbind(1:nrow(neighbors),neighbors)
  smooth_val = sapply(1:length(val),function(i){
    weighted.mean(val[neighbors[i,]],w=exp(-rho*(sqrt(apply((xgrid[neighbors[i,],]-matrix(xgrid[i,],
                                                                                          byrow=TRUE,
                                                                                          nrow=length(neighbors[i,]),
                                                                                          ncol=2))^2,1,sum))))*(abs(val[neighbors[i,]]-val[i])<cor_diff),
                  na.rm = TRUE)
  })
  return(smooth_val)
}

hard.threshold = function(val,thres){
  return(ifelse(val>thres,val,0))
}

smooth.multi.imgs = function(Y,xgrid,neighbors=NULL,rad_voxel=3,n_cor=0.9){
  sY = matrix(NA,nrow=nrow(Y),ncol=ncol(Y))
  for(i in 1:ncol(dat$Y_1)){
    sY[,i]=spatial.kernel.smooth.2D.eq.sp(Y[,i],xgrid,neighbors=neighbors,rad_voxel=rad_voxel,n_cor=n_cor,cor_diff = Inf)
  }
  return(sY)
}

find.spatial.group = function(cluster,xgrid,adj_dist=1,size=100){
  cluster_idx = which(cluster)
  adj_mat = (as.matrix(dist(xgrid[cluster_idx,]))<=adj_dist)
  group_idx = spatial.cluster(adj_mat)
  tab = table(group_idx)
  tab = tab[which(tab>size)]
  uni_group = as.numeric(names(tab))
  group = list()
  if(length(tab)>0){
  for(i in 1:length(tab)){
    group[[i]] = cluster_idx[group_idx == uni_group[i]]
  }
  }
  return(group)
}

create.power.two.regions.2D = function(xgrid,level=1){ 
  parent_region = create.four.regions.2D(xgrid,1)
  if(level>1){
    for(i in 1:(level-1)){
      uni_parent_region = unique(parent_region)
      child_region = rep(0,length=nrow(xgrid))
      for(j in 1:length(uni_parent_region)){
          temp_idx = parent_region==uni_parent_region[j]
          child_region[temp_idx] = create.four.regions.2D(xgrid[temp_idx,],1+(j-1)*4)
      }
      parent_region = child_region
    }
  }
  return(parent_region)
}


spatial.cluster = function(adj_mat){
  tmp = which(adj_mat, arr.ind=TRUE)
  edge_index_pair = tmp[tmp[,2]>tmp[,1],]
  edge = as.vector(t(edge_index_pair))
  gph = graph(edge, n=as.numeric(dim(adj_mat)[1]), directed=FALSE)
  max_connect = clusters(gph)$membership
  return(max_connect)
}

comp.cluster.mean = function(cluster,rho){
  G = max(cluster)
  temp_rho = rep(0,length=length(rho))
  for(g in 1:G)
    temp_rho = ifelse(cluster==g,mean(rho[cluster==g]),temp_rho)
  return(temp_rho)
}

standardize.data = function(X){
  newX = sapply(1:nrow(X), function(i) return((X[i,]-mean(X[i,]))/sd(X[i,])))
  return(newX)
}

compute.log.likelihood = function(X,Y,R_x_est,R_y_est,rho,modify=1e-10){
  n = nrow(X)
  p = ncol(X)
  Sigma = matrix(0,nrow=2*p,ncol=2*p)
  Sigma[1:p,1:p] = R_x_est
  Sigma[p+1:p,p+1:p] = R_y_est
  Sigma[1:p,p+1:p] = diag(p)*rho
  Sigma[p+1:p,1:p] = diag(p)*rho
  eigres = eigen(Sigma)
  if(min(eigres$values)<=0)
    loglike = -Inf
  else{
    XY = cbind(X,Y)
    temp = diag(1/sqrt(eigres$values))%*%t(eigres$vectors)%*%t(XY)
    loglike = 0.5*n*sum(log(eigres$values+modify))-0.5*sum(diag(t(temp)%*%temp))
  }
  return(loglike)
}

region.test = function(group_idx,dat_1,dat_2,rho_list = seq(-0.8,0.8,by=0.05)){
  X = dat_1[group_idx,]
  dim(X) = c(length(group_idx),ncol(dat_1))
  Y = dat_2[group_idx,]
  dim(Y) = c(length(group_idx),ncol(dat_2))
  X = standardize.data(X)
  Y = standardize.data(Y)
  p  = length(group_idx)
  R_x_est = diag(p)
  R_y_est = diag(p)
  loglike=sapply(rho_list, function(rho) compute.log.likelihood(X,Y,R_x_est,R_y_est,rho))
  Lambda = 2*(loglike[which.max(loglike)] - loglike[which(rho_list==0)])
  pvalue = pchisq(Lambda,df=1,lower.tail = FALSE)
  return(list(res=c(pvalue=pvalue,rho_est=rho_list[which.max(loglike)],Lambda=Lambda),
              loglike = cbind(rho_list=rho_list,loglike)))
}

sapply.pb <- function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}

lapply.pb <- function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}



region.cor.map.two.2D = function(Y_1,Y_2,xgrid,region){
  region_res = region.selection(Y_1,Y_2,xgrid,region)
  region_pos_idx = which(region_res$adj_region_pvalue<0.05 & region_res$region_value>0)
  region_neg_idx = which(region_res$adj_region_pvalue<0.05 & region_res$region_value<0)
  region_selection = rep(0,length=nrow(xgrid))
  region_selection[which(!is.na(match(region,region_res$uni_region[region_pos_idx])))] = 1
  region_selection[which(!is.na(match(region,region_res$uni_region[region_neg_idx])))] = -1
  cor_val = rep(0,length=nrow(xgrid))
  pvalue = rep(0,length=nrow(xgrid))
  adj_pvalue = rep(0,length=nrow(xgrid))
  for(i in 1:length(region_res$uni_region)){
    temp_idx = which(region==region_res$uni_region[i])
    cor_val[temp_idx] = region_res$region_value[i]
    pvalue[temp_idx] = region_res$region_pvalue[i]
    adj_pvalue[temp_idx] = region_res$adj_region_pvalue[i]
    
    
  }
  return(list(cor_type=region_selection,
              cor_val=cor_val,
              pvalue=pvalue,
              adj_pvalue = adj_pvalue))
}


voxel.cor.map.two.2D = function(Y_1,Y_2){
  voxel_cor = voxel.cor.est(Y_1,Y_2)
  voxel_cor_pvalue = voxel.cor.test(Y_1,Y_2)
  adj_voxel_pvalue = p.adjust(voxel_cor_pvalue)
  voxel_selection = ifelse((adj_voxel_pvalue<0.05) & (voxel_cor>0),1,0)
  voxel_selection = ifelse((adj_voxel_pvalue<0.05) & (voxel_cor<0),-1,voxel_selection)
  return(list(cor_type=voxel_selection,cor_val = voxel_cor,pvalue = voxel_cor_pvalue,
              adj_pvalue = adj_voxel_pvalue))
}

part.cor.map.two.2D = function(dat_1,dat_2,xgrid,rho_list = seq(-0.5,0.5,by=0.01),
                               adj_dist=1,size=100,thres_prob=0.85){
  
  voxel_cor = voxel.cor.est(dat_1,dat_2)
  
  smooth_cor = spatial.kernel.smooth.2D.eq.sp(voxel_cor,xgrid)
  thres_pos_cor = hard.threshold(smooth_cor, quantile(smooth_cor,prob=thres_prob,na.rm=TRUE))
  thres_neg_cor = -hard.threshold(-smooth_cor, quantile(-smooth_cor,prob=thres_prob,na.rm=TRUE))
  pos_cluster = thres_pos_cor>0
  neg_cluster = thres_neg_cor<0
  pos_group = find.spatial.group(pos_cluster,xgrid,adj_dist,size)
  neg_group = find.spatial.group(neg_cluster,xgrid,adj_dist,size)
  all_group = c(neg_group,pos_group)
  
  simple_est_all = sapply(1:length(all_group),function(i) mean(voxel_cor[all_group[[i]]]))
  mle_all = lapply.pb(1:length(all_group),
                     function(i) region.test(group_idx = all_group[[i]],dat_1,dat_2,rho_list))
  
  return(list(cor_val=voxel_cor,smooth_cor_est=smooth_cor,pos_cluster=pos_cluster,
              neg_cluster=neg_cluster,
              neg_group = neg_group,pos_group=pos_group,all_group=all_group,
              simple_est_all=simple_est_all,mle_all=mle_all))
}


spat.adapt.cor.map.two.2D = function(dat_1,dat_2,xgrid,rho_list = seq(-0.5,0.5,by=0.01),
                               adj_dist=1,size=100,thres_prob=0.85,neighbors=NULL,cor_diff=0.1){
  
  voxel_cor = voxel.cor.est(dat_1,dat_2)
  
  smooth_cor = spatial.kernel.smooth.2D.eq.sp(voxel_cor,xgrid,neighbors = neighbors,cor_diff=cor_diff)
  thres_pos_cor = hard.threshold(smooth_cor, quantile(smooth_cor,prob=thres_prob,na.rm=TRUE))
  thres_neg_cor = -hard.threshold(-smooth_cor, quantile(-smooth_cor,prob=thres_prob,na.rm=TRUE))
  pos_cluster = thres_pos_cor>0
  neg_cluster = thres_neg_cor<0
  pos_group = find.spatial.group(pos_cluster,xgrid,adj_dist,size)
  neg_group = find.spatial.group(neg_cluster,xgrid,adj_dist,size)
  selection = rep(0,length=nrow(xgrid))
  if(length(pos_group)>0){
  for(i in 1:length(pos_group))
    selection[pos_group[[i]]]=1
  }
  if(length(neg_group)>0){
  for(i in 1:length(neg_group))
    selection[neg_group[[i]]]=-1
  }
  
  
  return(list(cor_val=voxel_cor,smooth_cor_est=smooth_cor,pos_cluster=pos_cluster,
              neg_cluster=neg_cluster,
              neg_group = neg_group,pos_group=pos_group,cor_type=selection))
}

create.four.regions.2D = function(xgrid,start_idx = 1){
  center = apply(xgrid,2,mean)
  region = ifelse(xgrid[,1]<=center[1] & xgrid[,2]<=center[2],start_idx,0)
  region = ifelse(xgrid[,1]<=center[1] & xgrid[,2]>center[2],start_idx+1,region)
  region = ifelse(xgrid[,1]>center[1] & xgrid[,2]<=center[2],start_idx+2,region)
  region = ifelse(xgrid[,1]>center[1] & xgrid[,2]>center[2],start_idx+3,region)
  return(region)
}


region.selection = function(Y_1,Y_2,xgrid,region){
  uni_region = unique(region)
  region_value = rep(NA,length=length(uni_region))
  region_pvalue = rep(NA,length=length(uni_region))
  for(i in 1:length(uni_region)){
    temp = cor.test(apply(Y_1[region==uni_region[i],],2,mean),apply(Y_2[region==uni_region[i],],2,mean))
    region_value[i] = temp$estimate
    region_pvalue[i] = temp$p.value
  }
  adj_region_pvalue = p.adjust(region_pvalue)
  return(list(uni_region = uni_region, region_value=region_value,region_pvalue=region_pvalue,adj_region_pvalue=adj_region_pvalue))
}

create.res.table.2D = function(res){
  all_group = res$all_group
  mle_all = res$mle_all
  simple_est_all = res$simple_est_all
  all_tab = NULL
  
  for(i in 1:length(all_group)){
    all_tab = rbind(all_tab,mle_all[[i]]$res)
  }
  all_tab = cbind(all_tab,simple_est_all)
  
  print_p_value = function(p_value){
    return(ifelse(p_value<1e-8,"<1e-8",signif(p_value,digits=2)))
  }
  cluster_size = sapply(1:length(all_group),function(i) length(all_group[[i]]))
  print_tab = cbind(cluster_size,round(all_tab[,4],digit=3),round(all_tab[,2],digit=3),
                    print_p_value(p.adjust(all_tab[,1],method="BY")))
  colnames(print_tab) = c("cluster size","mean voxel cor","region cor","p value")
  idx = which(all_tab[,2]*all_tab[,4]>0)
  return(print_tab[idx,])
}

comp.region.center = function(xgrid,region){
  uni_region = unique(region)
  centers  = matrix(NA,nrow=length(uni_region),ncol=ncol(xgrid))
  for(i in 1:length(uni_region)){
    centers[i,] = apply(xgrid[region==uni_region[i],],2,mean)
  }
  return(cbind(uni_region,centers))
}

classify.accuracy = function(est_type,true_type,levels=NULL){
  if(is.null(levels))
    levels = unique(true_type)
  tab = table(factor(est_type,levels=levels), factor(true_type,levels=levels))
  TP = tab[1,1]
  FP = tab[1,2]
  TN = tab[2,2]
  FN = tab[2,1]
  Sensitivity = TP/(TP+FN)
  Specificity = TN/(TN+FP)
  PPV = TP/(TP+FP)
  NPV = TN/(TN+FN)
  Accuracy = (TP+TN)/sum(tab)
  FDR = 1 - PPV
  return(c(TP=TP,FP=FP,TN=TN,FN=FN,
           Sensitivity=Sensitivity,
           Specificity=Specificity,
           PPV=PPV,NPV=NPV,FDR=FDR,
           Accuracy=Accuracy))
}

print.mean.sd.tab = function(mean_tab,sd_tab){
  m = ncol(mean_tab)
  print_tab = NULL
  for(i in 1:m){
    print_tab = cbind(print_tab,sprintf("%.3f (%.3f)",mean_tab[,i],sd_tab[,i]))
  }
  
  colnames(print_tab) = colnames(mean_tab)
  return(print_tab)
}
