SBM_EMV = function(x,r,rho,mu,sigma2_1){
  
  P = dim(r)[2]
  N = dim(r)[1]
  
  eps = 10^(-3)
  max_iteration = 100
  Free_E0 = -Inf
  for(ite in 1:max_iteration){
    
    # r estimation
    lpsi = function(i,k){
      log_psi_out = 0
      for(j in 1:N){
        for(l in 1:P){
          if(i != j){
            log_psi_out = log_psi_out + r[j,l]*(-((x[i,j] - mu[k,l])^2)/(2*sigma2_1) - log(sqrt(2*pi*sigma2_1)))
          }
        }
      }
      return(log_psi_out)
    }
    
    log_psi = matrix(0,N,P)
    for(i in 1:N){
      for(k in 1:P){
        log_psi[i,k] = lpsi(i,k)
      }
    }
    log_psi0 = matrix(0,N,P)
    for(i in 1:N){
      log_psi0[i,] = log_psi[i,] - max(log_psi[i,])
    }
    r0 = matrix(0,N,P)
    for(i in 1:N){
      r_den = 0
      for(k in 1:P){
        r0[i,k] = rho[k]*exp(log_psi0[i,k])
        r_den = r_den + rho[k]*exp(log_psi0[i,k])
      }
      r0[i,] = r0[i,]/r_den
    }
    r = r0
    
    # rho estimation
    rho = as.matrix(apply(r,2,mean))
    
    # mu and sigma2 estimation
    mu_num = matrix(0,P,P)
    mu_den = matrix(0,P,P)
    for(k in 1:P){
      for(l in 1:P){
        for(i in 1:N){
          for(j in min((i+1),N):N){
            mu_num[k,l] = mu_num[k,l] + r[i,k]*r[j,l]*x[i,j]
            mu_den[k,l] = mu_den[k,l] + r[i,k]*r[j,l]
          }
        }
      }
    }
    mu = mu_num/mu_den
    
    sigma2_1_num = 0
    sigma2_1_den = 0
    for(k in 1:P){
      for(l in 1:P){
        for(i in 1:N){
          for(j in min((i+1),N):N){
            sigma2_1_num = sigma2_1_num + r[i,k]*r[j,l]*(x[i,j]-mu[k,l])^2
            sigma2_1_den = sigma2_1_den + r[i,k]*r[j,l]
          }
        }
      }
    }
    sigma2_1 = sigma2_1_num/sigma2_1_den
    
    # free energy
    Free_E = 0
    for(i in 1:N){
      for(k in 1:P){
        for(j in min((i+1),N):N){
          for(l in 1:P){
            Free_E = Free_E + r[i,k]*r[j,l]*(-((x[i,j] - mu[k,l])^2)/(2*sigma2_1) - 0.5*log(2*pi*sigma2_1))
          }
        }
        Free_E = Free_E + log(rho[k]^r[i,k]) - log(r[i,k]^r[i,k])
      }
    }
    
    if(abs((Free_E-Free_E0)/Free_E)<eps){
      break
    }
      Free_E0 = Free_E
      print(Free_E)
  }

  return(list("r"=r,"rho"=rho,"mu"=mu,"sig1"=sigma2_1,"Fe"=Free_E))
}
