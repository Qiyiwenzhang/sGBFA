library(MASS)
library(statmod)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(MCMCpack)
library(nloptr)
library(glmnet)
library(BayesLogit)
library(parallel)


##################### working graph ###########################
## x : 1 refers  working graph = true graph
##     2 refers  working graph = adding edges between pathways 
##     3 refers  working graph = random graph
##     4 refers  working graph = removing edges from the true graph 
working_graph <- function(x,pathway_list,H,data_dim,ind_s,ind_e,p) {
  graph = matrix(0,p,p)
  if(x==0){
    graph = matrix(0,p,p)
    return(graph)
  }else if(x==1){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]= rbinom(length(c_ind),1,0.7)
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = 0.5 * (graph + t(graph))
    return(graph)  
  }else if(x==2){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = 0.5 * (graph + t(graph))
    return(graph)  
  }else if(x==3){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
        graph_tmp[(node[nn]+1):sum(pathway_tmp[1:nn]),(node[nn]+1):sum(pathway_tmp[1:nn])] = rbinom((pathway_tmp[nn]-1)^2,1,0.3)
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }else if(x==4){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
        # within pathway 
        graph_tmp[(node[nn]+1):sum(pathway_tmp[1:nn]),(node[nn]+1):sum(pathway_tmp[1:nn])] = rbinom((pathway_tmp[nn]-1)^2,1,0.3)
        # across pathway 
        if(sum(pathway_tmp[1:nn])+1 < sum(pathway_tmp)){
          graph_tmp[node[nn]:sum(pathway_tmp[1:nn]),(sum(pathway_tmp[1:nn])+1):sum(pathway_tmp)] = rbinom(pathway_tmp[nn]*(sum(pathway_tmp)-sum(pathway_tmp[1:nn])),1,0.1) 
        }
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }
}



###########################  empty chain ##########################
empty_chain <- function(n,p,T,L,H,rho_ini,phi_ini,w_ini,z_ini,alpha_ini,tau_ini){
  chain_m   = matrix(0,ncol = p, nrow = T+2)
  
  chain_rho   = matrix(0,ncol = p*n, nrow = T+2)
  chain_rho[1,]=as.vector(t(rho_ini))
  
  chain_w   = matrix(0,ncol = p*L, nrow = T+2)
  chain_w[1,]=as.vector(t(w_ini))
  
  chain_phi   = matrix(0,ncol = H*L, nrow = T+2)
  chain_phi[1,]=as.vector(t(phi_ini))
  
  chain_tau   = matrix(0,ncol = p*L, nrow = T+2)
  chain_tau[1,]=as.vector(t(tau_ini))
  
  chain_z   = matrix(0,ncol = L*n, nrow = T+2)
  chain_z[1,]=as.vector(t(z_ini))
  
  chain_alpha =matrix(0,ncol = p*L, nrow = T+2)
  chain_alpha[1,]=as.vector(t(alpha_ini))
  
  list(chain_tau=chain_tau,chain_alpha=chain_alpha,chain_phi=chain_phi,chain_rho=chain_rho,chain_m=chain_m,chain_w=chain_w,chain_z=chain_z)
}



## propto density of alpha_l
# alpha_l : candidate value
# alpha_h : current value
# alpha   : p by 1 alpha vector 
# omega   : the precison matrix 
# tau_l   : tau[j,l]^2
# j : the j-th entry of alpha 


p_density_mmh <- function(alpha_l,alpha_h,tau_l,phi_l,w_l,alpha,omega,j,nu_1,nu_2,p) {
  # alpha_l : candiate
  # alpha_h : old 
  lambda_l = exp(alpha_l)
  lambda_h = exp(alpha_h)
  alp_l_aug = replace(alpha,j,alpha_l)
  omega_j = omega[j,]
  ans = (lambda_l/lambda_h)^(3/2)*exp(-(lambda_l-lambda_h)*(tau_l*phi_l+(w_l^2)*phi_l/tau_l)/2)*exp((-1/(2*nu_2))*(alpha_l-alpha_h)*omega_j%*%(alpha+alp_l_aug-2*nu_1*rep(1,p)))
  return(ans)
}


sample_quantile = function(x){
  probs = c(0.025,0.975)
  ans = quantile(x,probs = probs)
  return(ans)
} 


#### MCMC_algorithm #### 

## arguments of the function

# dt : {0,1,2,3} = {gaussian, binary, binomial,mix}
# T : the number of iterations of MCMC
# L : the number of factors  
# p : the number of features
# n : the number of subjects 
# X : data 


# nu_1  :  mean 
# nu_2  :  variance 
# Sigma :  variance for each entry of m
# Q : covariance matrix for the proposal density
# eta : hyper parameter specifying the prior of omega
# eps : hyper parameter specifying the prior of omega 

# rho_temp : initial value for rho 
# tau_temp : initial value for tau 
# omega_temp : initial value for omega
# inv_omega_temp : initial value for inverse omega
# w_temp : initial value for w
# z_temp : initial value for z
# alpha_temp : initial value for 
# m_temp : initial value for

# th  : threshold to determine the degree of freedom (of W)

##########################################

## values of the function
# w_est : estimation of W
# z_est : estimation of Z
# m_est : estimation of m 
# df    : degree of freedom 
# BIC   : BIC criterion 
# chain_alpha : markov chain of alpha 
# chain_m     : markov chain of m 
# chain_w     : markov chain of W
# chain_z     : markov chain of Z 


BFGA_MCMC <- function(dt,outcome,T,L,p,n,X,trials,ind_s,ind_e,data_dim,a_lambda,b_lambda,nu_1,nu_2,Sigma,Q,eta,eps,a_phi,b_phi,phi_temp,rho_temp,tau_temp,omega_temp,inv_omega_temp,w_temp,z_temp,alpha_temp,m_temp,start,end,H,mu,chain_rho,chain_alpha,chain_tau,chain_w,chain_z,chain_m,chain_phi){
  
  
  
  # preparation based on the data type 
  
  no_dt = unique(dt)
  dt_gau = which(dt==0)
  dt_ber = which(dt==1)
  dt_bin = which(dt==2)
  block_id <- rep(seq_len(H), times = data_dim)
 
  # only depend on data types 
  if(length(no_dt)==1){
    if(dt[1]==0){
      psi = X
      ## phi: prior parameter for rho of gau  
      phi =  as.matrix(rep(1,p))
      kappa = matrix(0,nrow = p,ncol = n)
      rowsum_kappa = rep(0,p)
    }else if(dt[1]==1){
      psi = matrix(0,nrow = p,ncol = n)
      b = matrix(1,ncol = n,nrow = p) # binary data
      kappa = X-b/2
      rowsum_kappa= rowSums(kappa)
    }else if(dt[1]==2){
      psi = matrix(0,nrow = p,ncol = n)
      b = trials
      kappa = X-b/2
      rowsum_kappa= rowSums(kappa)
    }
  } else{
    psi =  matrix(0,nrow = p,ncol = n)
    psi[dt_gau,] = X[dt_gau,]
    b  = matrix(1,ncol = n,nrow = p)
    b[dt_bin,] = trials[dt_bin,]
    kappa = X - b/2
    kappa[dt_gau,] = 0
    rowsum_kappa= rowSums(kappa)
    ## phi: prior parameter for rho of gau  
    phi = as.matrix(rep(1,p))
  }
  
  
  for (t in 1:T) {
    rowsum_rho = rowSums(rho_temp)
    row_rho_psi = rowSums(rho_temp*psi)
    
    
    ###################################### update m ######################################
    for (j in 1:p) {
      covar_m = (rowsum_rho[j]+Sigma[j]^-1)^-1
      mean_m  = covar_m*(rowsum_kappa[j]+row_rho_psi[j]- (w_temp[j,]%*%z_temp)%*%rho_temp[j,])
      m_temp[j] = rnorm(1,mean_m,sqrt(covar_m))
    }
    
    chain_m[t+1,]=as.vector(m_temp)
    
    ####################################### update rho ####################################
    for (j in 1:p) {
      if(j %in% dt_gau){
        # rho_j
        r = phi[j,1]+sum((X[j,]-(m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])))^2)
        rho_j = rgamma(1,shape=(phi[j,1]+n)/2 ,rate = r/2 )*rep(1,n)
        rho_temp[j,] = rho_j  
      }else{
        mu_j = (m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])) 
        for (i in 1:n) {
          rho_temp[j,i] = rpg.devroye(1,b[j,i],mu_j[i])  
        }
      }
    }
    
    chain_rho[t+1,]=as.vector(t(rho_temp))
    
    ####################################### generate alpha_l ####################################### 
    
    ## set up outcome range: the last modality 
    if(outcome==1){
      outcome_j_range = seq(ind_s[H],ind_e[H],1)  
    }else{
      outcome_j_range = c(0)
    }
    
    if (length(no_dt)==1){
      for (l in 1:L) {
        for (j in 1:p) {
          if(j %in% outcome_j_range){
            lambda_shape = a_lambda + 3/2
            lambda_rate  = b_lambda + (tau_temp[j,l])/2 + (w_temp[j,l])^2/(2*tau_temp[j,l]) 
            alpha_temp[j,l] = log(rgamma(1,shape=lambda_shape,rate=lambda_rate)) # change lambda to log-scale for the outcome modality 
          }else{
            tmp_h <- block_id[j]
            candi = rnorm(1,alpha_temp[j,l],1)
            comp =  as.numeric(p_density_mmh(candi,alpha_temp[j,l],tau_temp[j,l],phi_temp[tmp_h,l],w_temp[j,l],alpha_temp[,l],omega_temp,j,nu_1,nu_2,p))   
            if(is.na(comp)){
              comp=0
            }
            # print(comp)
            prob = min(comp,1)
            u = runif(1,0,1)
            #  print(u)
            if(u<prob){
              alpha_temp[j,l]=candi
            }
          }
        }
      }
    }else{
      for (l in 1:L) {
        for (j in 1:p) {
          
          if(j %in% outcome_j_range){
            lambda_shape = a_lambda + 3/2
            lambda_rate  = b_lambda + (tau_temp[j,l])/2 + (w_temp[j,l])^2/(2*tau_temp[j,l]) 
            alpha_temp[j,l] = log(rgamma(1,shape=lambda_shape,rate=lambda_rate)) # change lambda to log-scale for the outcome modality 
          }else{
          ## assign view-specific 
          if(j %in% dt_gau){
            nu_1.viewsep = nu_1[1]
            nu_2.viewsep = nu_2[1]
          }else if (j %in% dt_ber){
            nu_1.viewsep = nu_1[2]
            nu_2.viewsep = nu_2[2]
          }else if(j %in% dt_bin){
            nu_1.viewsep = nu_1[3]
            nu_2.viewsep = nu_2[3]
          }
          
          tmp_h <- block_id[j]
          candi = rnorm(1,alpha_temp[j,l],1)
          comp =  as.numeric(p_density_mmh(candi,alpha_temp[j,l],tau_temp[j,l],phi_temp[tmp_h,l],w_temp[j,l],alpha_temp[,l],omega_temp,j,nu_1.viewsep,nu_2.viewsep,p))   
          if(is.na(comp)){
            comp=0
          }
          # print(comp)
          prob = min(comp,1)
          u = runif(1,0,1)
          #  print(u)
          if(u<prob){
            alpha_temp[j,l]=candi
          }
          }
        }
      }
    }
    
    chain_alpha[t+1,]=as.vector(t(alpha_temp))
    
    
    ## calculate A  
    
    if (length(no_dt)==1){
      A = eta*(diag(eps,p) + 1) + (alpha_temp-nu_1)%*%t(alpha_temp-nu_1)/nu_2
      
    } else{
      nu_1.vec = rep(nu_1,each=p/3) 
      nu_2.vec = rep(nu_2,each=p/3)
      A = eta*(diag(eps,p) + 1) + ((alpha_temp-nu_1.vec)/sqrt(nu_2.vec)) %*% t(((alpha_temp-nu_1.vec)/sqrt(nu_2.vec)))
    }
    
    ####################################### update Omega   #######################################
    
    for (j in 1:p) {
      ## j-th diagonal element and j-th column excluding the diag 
      A_jj = A[j,j]     
      a_j = A[-j,j]   
      
      # construct the inverse of Omega_11 via inverse of Omega 
      inv_omega_11 = inv_omega_temp[-j,-j]
      inv_omega_12 = inv_omega_temp[-j,j]
      inv_omega_22 = inv_omega_temp[j,j]
      
      O_11 = inv_omega_11 - inv_omega_12%*%t(inv_omega_12)/inv_omega_22
      
      ## find the index of non-zeros of j-th col 
      nz_ind = which(omega_temp[,j]!=0) # include diagonal
      
      ## remove the diag index: index==j
      nz_ind = nz_ind[nz_ind!=j]        # exclude diagonal
      
      
      ## see if there is any non-zero entry except for the diag
      ## if length(nz_ind)=0: only generate the diag
      if(length(nz_ind)!=0){
        
        ## define the block matrices (mean part)
        ## i) select nz_ind rows and remove the j-th col 
        ## ii) remove the (nz_ind,j) rows and the j-th col 
        ## construct the precision matrix for \omega_12^{(1)}
        
        if(length(nz_ind)>1){
          prec_nz = inv_omega_temp[nz_ind, nz_ind]-inv_omega_temp[nz_ind,j]%*%t(inv_omega_temp[nz_ind,j])/inv_omega_temp[j,j]
          prec_nz_sqrt = chol(prec_nz)  # upper tri matrix 
          mid = rnorm(length(nz_ind),0,1/sqrt(A_jj))
          pt_mean_o1 = (-1/A_jj)*A[nz_ind,j]
          pt_mean_o2 = forwardsolve(t(prec_nz_sqrt),pt_mean_o1)
          mean_nz = backsolve(prec_nz_sqrt,pt_mean_o2)
          non_zero_omega = backsolve(prec_nz_sqrt,mid) + mean_nz 
          
          ## generate the diag entry
          
          xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
          partial = prec_nz_sqrt%*%non_zero_omega
          diag_omega = xi + t(partial)%*%partial
          
        }else if(length(nz_ind)==1){
          prec_nz = inv_omega_temp[nz_ind, nz_ind]-inv_omega_temp[nz_ind,j]%*%t(inv_omega_temp[nz_ind,j])/inv_omega_temp[j,j]
          mean_nz = (-1/(A_jj*prec_nz))*A[nz_ind,j]
          non_zero_omega = rnorm(1,mean_nz, 1/sqrt(prec_nz*A_jj))
          
          ## generate the diag entry
          
          xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
          diag_omega = xi + (non_zero_omega^2)*prec_nz
        }
        
        ## fill in the non-zero entries with the new sample
        ## both col and row d
        omega_temp[nz_ind,j] = non_zero_omega   # col 
        omega_temp[j,nz_ind] = non_zero_omega   # row
        omega_temp[j,j]  = diag_omega
      }else {
        
        ## when each entry = 0 except for the diag  
        ## we only generate the xi with Gamma
        xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
        diag_omega = xi
        omega_temp[j,j]  = diag_omega
      }
      
      # update inv_omega_temp 
      inv_omega_temp[j,j] = 1/xi
      pp = O_11%*%omega_temp[-j,j]
      inv_omega_11 = O_11 + inv_omega_temp[j,j]*pp%*%t(pp)
      inv_omega_12 = -inv_omega_11%*%omega_temp[-j,j]/omega_temp[j,j]
      
      inv_omega_temp[-j,-j] = inv_omega_11
      inv_omega_temp[-j,j] = inv_omega_12
      inv_omega_temp[j,-j] = inv_omega_12 
      
    }
    
    ####################################### generate z_i   #######################################
    
    for (i in 1:n) {
      pt_cov = t(w_temp)%*%(w_temp*rho_temp[,i])
      diag(pt_cov) = diag(pt_cov)+ rep(1,L)
      chol_pt_cov = chol(pt_cov) # upper matrix 
      pt_mean_z1 = t(w_temp)%*%(rho_temp[,i]*(psi[,i]-m_temp)+kappa[,i])
      pt_mean_z2 = forwardsolve(t(chol_pt_cov),pt_mean_z1)
      mean_z_i  = backsolve(chol_pt_cov,pt_mean_z2)
      z_temp[,i] = backsolve(chol_pt_cov,rnorm(L,0,1)) + mean_z_i
    }
    chain_z[t+1,]=as.vector(t(z_temp))
    
    ############################################# update phi   #######################################
    
    if(outcome==1){
      for (h in 1:(H-1)) {
        ind_h_s = ind_s[h]
        ind_h_e = ind_e[h]
        p_h = data_dim[h]
        for (l in 1:L) {
          phi_shape = a_phi + 3*p_h/2
          phi_rate  = b_phi + sum((tau_temp[ind_h_s:ind_h_e,l])*exp(alpha_temp[ind_h_s:ind_h_e,l]))/2+
            sum(exp(alpha_temp[ind_h_s:ind_h_e,l])*(w_temp[ind_h_s:ind_h_e,l])^2/(2*(tau_temp[ind_h_s:ind_h_e,l])))
          phi_temp[h,l] = rgamma(1,shape=phi_shape, rate=phi_rate)
        }
      }
      
    }else{
      for (h in 1:H) {
        ind_h_s = ind_s[h]
        ind_h_e = ind_e[h]
        p_h = data_dim[h]
        for (l in 1:L) {
          phi_shape = a_phi + 3*p_h/2
          phi_rate  = b_phi + sum((tau_temp[ind_h_s:ind_h_e,l])*exp(alpha_temp[ind_h_s:ind_h_e,l]))/2+
            sum(exp(alpha_temp[ind_h_s:ind_h_e,l])*(w_temp[ind_h_s:ind_h_e,l])^2/(2*(tau_temp[ind_h_s:ind_h_e,l])))
          phi_temp[h,l] = rgamma(1,shape=phi_shape, rate=phi_rate)
        }
      }
    }
    
    chain_phi[t+1,]=as.vector(t(phi_temp))
    
    
    # update w_j (include regression coefficients) & tau
    #################################################################
    for (j in 1:p) {
      
      tmp_h <- block_id[j]
      
      cov_pt = z_temp%*%(t(z_temp)*rho_temp[j,])
      diag(cov_pt)= diag(cov_pt)+ (1/tau_temp[j,])*(exp(alpha_temp[j,])*phi_temp[tmp_h,])
      chol_cov_pt = chol(cov_pt)
      
      pt_mean_w1 = t(t(z_temp)*rho_temp[j,])%*%(psi[j,]-m_temp[j]*rep(1,n)+(rho_temp[j,]^-1)*kappa[j,])
      pt_mean_w2 = forwardsolve(t(chol_cov_pt),pt_mean_w1)
      mean_w_j = backsolve(chol_cov_pt,pt_mean_w2)
      w_temp[j,] = backsolve(chol_cov_pt,rnorm(L,0,1)) + mean_w_j
      tau_temp[j,]=1/rinvgauss(L,1/abs(w_temp[j,]),exp(alpha_temp[j,])*phi_temp[tmp_h,])
    }
    chain_w[t+1,]=as.vector(t(w_temp))
    chain_tau[t+1,]=as.vector(t(tau_temp))
    #print(paste0("Task Progress: ", t/T ))
  }
  
  
  ################################## Result Summary ################################## 
  likelihood_chain = c() 
  mu_est_t=matrix(0,nrow=p,ncol=n)
  
  for (t in start:end) {
    prec_mat_t = matrix(chain_rho[t,],nrow=p,ncol=n,byrow=T)
    w_est_t = matrix(chain_w[t,],nrow=p,ncol=L,byrow=TRUE)
    z_est_t = matrix(chain_z[t,],nrow=L ,ncol=n,byrow=TRUE)
    m_est_t = matrix(chain_m[t,],nrow=p ,ncol=1,byrow=T)
    
    mu_est_t = mu_est_t + w_est_t%*%z_est_t + as.vector(m_est_t)
    mu_est_t.current = w_est_t%*%z_est_t + as.vector(m_est_t)
    
    # matrix form 
    likelihood_gau_t  = 0.5*log(prec_mat_t[dt_gau,])-0.5*log(2*pi)-0.5*prec_mat_t[dt_gau,]*(X[dt_gau,]-mu_est_t.current[dt_gau,])^2
    likelihood_ber_t = log(choose(trials[dt_ber,],X[dt_ber,]))+mu_est_t.current[dt_ber,]*X[dt_ber,]-log(1+exp(mu_est_t.current[dt_ber,]))*trials[dt_ber,]   
    likelihood_bin_t = log(choose(trials[dt_bin,],X[dt_bin,]))+mu_est_t.current[dt_bin,]*X[dt_bin,]-log(1+exp(mu_est_t.current[dt_bin,]))*trials[dt_bin,]   
    # chain of sum of likelihood 
    likelihood_chain = c(likelihood_chain,sum(likelihood_gau_t)+sum(likelihood_ber_t)+sum(likelihood_bin_t))
  }
  
  ################################## calculate RE and prediction error (in sample) ################################## 
  mu_est =mu_est_t/(end-start+1)
  if(outcome==1){
    # dim reduction error 
    mu_est_mod = mu_est[-outcome_j_range,]
    mu_mod  = mu[-outcome_j_range,]
    error_mod = norm(mu_est_mod-mu_mod,type = "F")/norm(mu_mod,type = "F")  
    # outcome prediction (in-sample)
    mu_est_outcome = mu_est[outcome_j_range,]
    X_outcome = X[outcome_j_range,]
    error_outcome = mean((mu_est_outcome-X_outcome)^2)
    error = c(error_mod, error_outcome)
  }else{
    error = norm(mu_est-mu,type = "F")/norm(mu,type = "F")  
  }
  
  
  # viewspecific error 
  if (length(no_dt)!=1){
    error.viewspe = rep(0,3)
    error.viewspe[1] =  norm(mu_est[dt_gau,]-mu[dt_gau,],type = "F")/norm(mu[dt_gau,],type = "F")
    error.viewspe[2] =  norm(mu_est[dt_ber,]-mu[dt_ber,],type = "F")/norm(mu[dt_ber,],type = "F")
    error.viewspe[3] =  norm(mu_est[dt_bin,]-mu[dt_bin],type = "F")/norm(mu[dt_bin,],type = "F")
  }else{
    error.viewspe = rep(0,3)
  }
  
  
  ################################## Degree of Freedom ################################## 
  interval = apply(chain_w[start:end,],2,sample_quantile)
  df = sum(interval[1,]*interval[2,]>0)
  
  ################################## BIC and DIC ################################## 
  
  prec_mat = matrix(apply(chain_rho[start:end,], 2, mean),nrow=p,ncol=n,byrow=T)
  likelihood_gau  = 0.5*log(prec_mat[dt_gau,])-0.5*log(2*pi)-0.5*prec_mat[dt_gau,]*(X[dt_gau,]-mu_est[dt_gau,])^2
  likelihood_ber = log(choose(trials[dt_ber,],X[dt_ber,]))+mu_est[dt_ber,]*X[dt_ber,]-log(1+exp(mu_est[dt_ber,]))*trials[dt_ber,]   
  likelihood_bin = log(choose(trials[dt_bin,],X[dt_bin,]))+mu_est[dt_bin,]*X[dt_bin,]-log(1+exp(mu_est[dt_bin,]))*trials[dt_bin,]   
  
  likelihood_pointest =(sum(likelihood_gau)+sum(likelihood_ber)+sum(likelihood_bin))
  
  # mean of likelihood 
  mean_like = mean(likelihood_chain)  
  dic = -2*likelihood_pointest + 2*2*(likelihood_pointest-mean_like)
  bic = -2*likelihood_pointest + log(n)*df
  
  
  ## estimation
  w_est = apply(chain_w[start:end,],2,mean)
  z_est = apply(chain_z[start:end,],2,mean)
  m_est = apply(chain_m[start:end,],2,mean)
  
  ###################################################################  
  list(error=error,error.viewspe=error.viewspe,interval=interval,df=df,bic=bic,dic=dic,w_est=w_est,z_est=z_est,m_est = m_est, mu_est=mu_est)
}


## helper: run one parameter combo
run_one_combo <- function(outcome,L, nu_1, nu_2,a_lambda,b_lambda, idx, p, n, H, T,
                          X, mu, trials, dt_vec,
                          ind_s, ind_e, data_dim,
                          graph, Sigma, Q, eta, eps,
                          a_phi, b_phi, start, end) {
  tryCatch({
    # init values ...
    rho_ini <- matrix(1, nrow=p, ncol=n)
    tau_ini <- matrix(1, nrow=p, ncol=L)
    set.seed(1234+idx)
    w_ini <- matrix(rnorm(L*p,0,1),nrow=p,ncol=L)
    z_ini <- matrix(rnorm(L*n,0,1),nrow=L,ncol=n)
    
    omega_ini <- diag(1, nrow=p, ncol=p)
    ind <- which(matrix(as.logical(graph), nrow=dim(graph)[1]), arr.ind=TRUE)
    if (nrow(ind) != 0) {
      for (k in 1:nrow(ind)) {
        if (ind[k,1] > ind[k,2]) omega_ini[ind[k,1], ind[k,2]] <- 0.05
      }
    }
    omega_temp     <- as.matrix(forceSymmetric(omega_ini, uplo='L'))
    inv_omega_temp <- solve(omega_temp)
    
    alpha_ini <- matrix(nu_1, nrow=p, ncol=L)
    phi_ini   <- matrix(1, nrow=H, ncol=L)
    
    mcmc_box <- empty_chain(n,p,T,L=L,H=H,
                            rho_ini=rho_ini,phi_ini=phi_ini,
                            w_ini=w_ini,z_ini=z_ini,
                            alpha_ini=alpha_ini,tau_ini=tau_ini)
    
    chain_rho   <- mcmc_box$chain_rho
    chain_w     <- mcmc_box$chain_w
    chain_z     <- mcmc_box$chain_z
    chain_alpha <- mcmc_box$chain_alpha
    chain_m     <- mcmc_box$chain_m
    chain_tau   <- mcmc_box$chain_tau
    chain_phi   <- mcmc_box$chain_phi
    
    w_temp     <- matrix(as.numeric(chain_w[1,]),nrow=p,ncol=L,byrow=TRUE)
    z_temp     <- matrix(as.numeric(chain_z[1,]),nrow=L,ncol=n,byrow=TRUE)
    rho_temp   <- matrix(as.numeric(chain_rho[1,]),nrow=p,ncol=n,byrow=TRUE)
    alpha_temp <- matrix(as.numeric(chain_alpha[1,]),nrow=p,ncol=L,byrow=TRUE)
    phi_temp   <- matrix(as.numeric(chain_phi[1,]),nrow=H,ncol=L,byrow=TRUE)
    m_temp     <- rep(0,p)
    
    trial_result <- BFGA_MCMC(dt=dt_vec,outcome = outcome, T=T,L=L,p=p,n=n,
                              X=X,trials=trials,
                              ind_s=ind_s,ind_e=ind_e,data_dim=data_dim,
                              a_lambda=a_lambda,b_lambda=b_lambda,
                              nu_1=nu_1,nu_2=nu_2,
                              Sigma=Sigma,Q=Q,eta=eta,eps=eps,
                              a_phi=a_phi,b_phi=b_phi,phi_temp=phi_temp,
                              rho_temp=rho_temp,tau_temp=tau_temp,
                              omega_temp=omega_temp,inv_omega_temp=inv_omega_temp,
                              w_temp=w_temp,z_temp=z_temp,alpha_temp=alpha_temp,
                              m_temp=m_temp,
                              start=start,end=end,H=H,mu=mu,
                              chain_rho=chain_rho,chain_alpha=chain_alpha,
                              chain_tau=chain_tau,chain_w=chain_w,
                              chain_z=chain_z,chain_m=chain_m,chain_phi=chain_phi)
    
    trial_result$L <- L
    trial_result
  }, error=function(e) {
    message("Combo failed: L=",L," nu1=",nu_1," nu2=",nu_2," : ",e$message)
    NULL
  })
}

sGBFA_sim <- function(idx,outcome, T, start, end, g, path=getwd()) {
  sim_list <- vector("list", s)
  data_1 = save_data_mask[[idx]]
  W = data_1$V
  Z = t(data_1$U)
  X = data_1$X
  mu = W %*% Z
  outcome = outcome
  ## dim parameter 
  n = dim(X)[2]
  true_L = dim(W)[2]
  H=data_1$H
  # p1 = data_1$p1
  # p2 = data_1$p2
  # p3 = data_1$p3
  # p4 = data_1$p4
  # p5 = data_1$p5
  p = data_1$p
  dt_vec = data_1$dt_vec
  data_dim = data_1$data_dim
  ind_s = data_1$ind_s
  ind_e = data_1$ind_e
  trials = data_1$trial_train
  pathway_list= data_1$pathway_list
  
  ## dim parameter 
  
  for (sim in 1:s) {
    dt_vec = data_1$dt_vec
    
    # data_tbs = list('true_W'= W,'true_Z'=Z,'true_X'=X)
    #####################################################
    
    Sigma = rep(0.01,p)  # variance for each entry of m
    
    Q = diag(4,nrow = p,ncol = p) # proposal density
    
    ######## working graph####################################
    
    if(outcome==1){
      graph = working_graph(g,pathway_list=pathway_list,H=H-1,data_dim =data_dim,ind_s=ind_s,ind_e=ind_e,p=sum(data_dim[1:(H-1)])) 
    }else{
      graph = working_graph(g,pathway_list=pathway_list,H=H,data_dim =data_dim,ind_s=ind_s,ind_e=ind_e,p=p) 
    }
   
    # data_tbs[['G']] = graph 
    ########################################################
    ## fixed hyper-para 
    a_phi = 1
    b_phi = 1
    eta = 10
    eps = 0.2
    a_lambda = 1
    b_lambda = 1
    ## gaussian  data: tuning para 
    grid_L  = c(true_L-1,true_L,true_L+1)
    ## gaussian data
    grid_nu_1 = c(0,0.5,1)   # mean 
    grid_nu_2 = c(0.5) # variance 
    
    no_dt = unique(dt_vec)
    if(length(no_dt)==1){
      params <- expand.grid(
        l = seq_along(grid_L),
        j = seq_along(grid_nu_1),
        i = seq_along(grid_nu_2)
      )
    }else{
      params <- expand.grid(
        l = seq_along(grid_L),
        j = seq_along(grid_nu_1[1,]),
        i = seq_along(grid_nu_2[1,])
      )
    }
    
    # ## parallel loop across parameter combinations
    trial_list_full <- mclapply(1:nrow(params), function(ii) {
      l <- params$l[ii]
      j <- params$j[ii]
      i <- params$i[ii]
      
      ## extract hyper-parameters
      nu_1 <- if (is.matrix(grid_nu_1)) grid_nu_1[, j] else grid_nu_1[j]
      nu_2 <- if (is.matrix(grid_nu_2)) grid_nu_2[, i] else grid_nu_2[i]
      L    <- grid_L[l]
      
      message(paste0("Running combo: L=", L, ", nu1=", nu_1, ", nu2=", nu_2))
      
      
      ############################ Initialization ############################
      print(paste0("L: ", L))
      
      ########################################################
      rho_ini = matrix(1,nrow=p,ncol=n )
      tau_temp = matrix(1,nrow=p,ncol=L )
      tau_ini  = tau_temp
      
      # set seed for initial value (modifed)
      set.seed(1234+idx)
      w_ini = matrix(rnorm(L*p,0,1),nrow=p,ncol=L)
      z_ini = matrix(rnorm(L*n,0,1),nrow=L,ncol=n)
      
      ## omega should be positive definite
      omega_ini = diag(1,nrow=p,ncol=p)  # compatible with graph
      
      ind = which(matrix(as.logical(graph), nrow=dim(graph)[1]), arr.ind=TRUE)
      if(dim(ind)[1]!=0){
        for (i in 1:dim(ind)[1]) {
          if(ind[i,1]>ind[i,2]){
            r_ind = ind[i,1]
            c_ind = ind[i,2]
            omega_ini[r_ind,c_ind] = 0.05
          }
        }
      }
      omega_temp = as.matrix(forceSymmetric(omega_ini,uplo = 'L'))
      inv_omega_temp = solve(omega_temp)
      
      
      if(is.matrix(grid_nu_1)){
        alpha_ini = matrix(rep(nu_1,each=p/H),nrow=p,ncol=L )
      }else{
        alpha_ini = matrix(nu_1,nrow=p,ncol=L )
        
      }
      phi_ini = matrix(1,nrow=H,ncol=L)
      mcmc_box = empty_chain(n,p,T,L=L, H=H,rho_ini=rho_ini,phi_ini=phi_ini,w_ini=w_ini,z_ini=z_ini,alpha_ini=alpha_ini,tau_ini=tau_ini)
      chain_rho = mcmc_box$chain_rho
      chain_w = mcmc_box$chain_w
      chain_z = mcmc_box$chain_z
      chain_alpha = mcmc_box$chain_alpha
      chain_m = mcmc_box$chain_m
      chain_tau = mcmc_box$chain_tau
      chain_phi = mcmc_box$chain_phi
      w_temp   = matrix(as.numeric(chain_w[1,]),nrow=p ,ncol=L,byrow=T)
      z_temp   = matrix(as.numeric(chain_z[1,]),nrow=L ,ncol=n,byrow=T)
      rho_temp   = matrix(as.numeric(chain_rho[1,]),nrow=p ,ncol=n,byrow=T)
      alpha_temp = matrix(as.numeric(chain_alpha[1,]),nrow=p ,ncol=L,byrow=T)
      phi_temp = matrix(as.numeric(chain_phi[1,]),nrow=H ,ncol=L,byrow=T)
      m_temp = rep(0,p)
      
      trial_result = BFGA_MCMC(dt=dt_vec,outcome = outcome, T=T,L=L,p=p,n=n,X,trials=trials,ind_s=ind_s,ind_e=ind_e,data_dim=data_dim,a_lambda,b_lambda,nu_1=nu_1,nu_2=nu_2,
                               Sigma,Q,eta,eps,a_phi=a_phi,b_phi=b_phi,phi_temp,
                               rho_temp,tau_temp,omega_temp,inv_omega_temp,w_temp,z_temp,alpha_temp,m_temp,start=start,end=end,H=H,mu=mu
                               ,chain_rho=chain_rho,chain_alpha=chain_alpha,chain_tau=chain_tau,chain_w=chain_w,chain_z=chain_z,chain_m=chain_m,chain_phi=chain_phi)
      trial_result$L = L
      trial_result
      
    }, mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset=1)))
    ########################################################
    # str(trial_list_full)
    ########################################################
    lowest_dic <- Inf
    best_id <- NULL
    gridL <- c()
    gride <- c()
    griddic <- c()

    for (i in seq_along(trial_list_full)) {
      current_dic <- trial_list_full[[i]]$dic

      gridL <- c(gridL, trial_list_full[[i]]$L)
      gride <- c(gride, trial_list_full[[i]]$error)
      griddic <- c(griddic, trial_list_full[[i]]$dic)

      if (current_dic < lowest_dic) {
        lowest_dic <- current_dic
        best_id <- i
      }
    }
    if (lowest_dic == Inf) {
      sim_list[[sim]] <- trial_list_full[[1]]
      sim_list[[sim]]$gridL <- gridL
      sim_list[[sim]]$gride <- gride
      sim_list[[sim]]$griddic <- griddic
    }
    else{
      sim_list[[sim]] <- trial_list_full[[best_id]]
      sim_list[[sim]]$gridL <- gridL
      sim_list[[sim]]$gride <- gride
      sim_list[[sim]]$griddic <- griddic
    }
  }

  saveRDS(trial_list_full, file=paste0(path,'/grid/grid_list_',g,'_',idx,'.rds'))
  saveRDS(sim_list, file=paste0(path,'/newsim_list_',g,'_',idx,'.rds'))
}



