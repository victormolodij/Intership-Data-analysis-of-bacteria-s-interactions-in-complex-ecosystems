SLBM_EMV = function(x,y,r,t,rho,tau,mu,nu,sigma2_1,sigma2_2){

  eps = 10^(-3)
  max_iteration = 100
  Free_E0 = -Inf
  for(ite in 1:max_iteration){
    
    # r,t estimation
    lpsi = function(i,k){
      log_psi_out = 0
      for(j in 1:N){
        for(l in 1:P){
          if(i != j){
            log_psi_out = log_psi_out + r[j,l]*(-((x[i,j] - mu[k,l])^2)/(2*sigma2_1) - 0.5*log(2*pi*sigma2_1))
          }
        }
      }
      for(g in 1:M){
        for(h in 1:Q){
            log_psi_out = log_psi_out + t[g,h]*(-((y[i,g] - nu[k,h])^2)/(2*sigma2_2) - 0.5*log(2*pi*sigma2_2))
        }
      }
      #return(exp(log_psi_out))
      return(log_psi_out)
    }
    
    lphi = function(g,h){
      log_phi_out = 0
      for(i in 1:N){
        for(k in 1:P){
          log_phi_out = log_phi_out + r[i,k]*(-((y[i,g] - nu[k,h])^2)/(2*sigma2_2) - 0.5*log(2*pi*sigma2_2))
        }
      }
      #return(exp(log_phi_out))
      return(log_phi_out)
    }
    
    log_psi = matrix(0,N,P)
    for(i in 1:N){
      for(k in 1:P){
        log_psi[i,k] = lpsi(i,k)
      }
    }
    log_psi_max = max(log_psi)
    log_psi0 = log_psi - log_psi_max
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
    
    log_phi = matrix(0,M,Q)
    for(g in 1:M){
      for(h in 1:Q){
        log_phi[g,h] = lphi(g,h)
      }
    }
    log_phi_max = max(log_phi)
    log_phi0 = log_phi - log_phi_max
    t0 = matrix(0,M,Q)
    for(g in 1:M){
      t_den = 0
      for(h in 1:Q){
        t0[g,h] = tau[h]*exp(log_phi0[g,h])
        t_den = t_den + tau[h]*exp(log_phi0[g,h])
      }
      t0[g,] = t0[g,]/t_den
    }
    t = t0
    
    # rho,tau estimation
    rho = as.matrix(apply(r,2,mean))
    tau = as.matrix(apply(t,2,mean))
    
    
    # mu,nu,sigma2 estimation
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
    
    nu_num = matrix(0,P,Q)
    nu_den = matrix(0,P,Q)
    for(k in 1:P){
      for(h in 1:Q){
        for(i in 1:N){
          for(g in 1:M){
            nu_num[k,h] = nu_num[k,h] + r[i,k]*t[g,h]*y[i,g]
            nu_den[k,h] = nu_den[k,h] + r[i,k]*t[g,h]
          }
        }
      }
    }
    nu = nu_num/nu_den
    
    sigma2_2_num = 0
    sigma2_2_den = 0
    for(k in 1:P){
      for(h in 1:Q){
        for(i in 1:N){
          for(g in 1:M){
            sigma2_2_num = sigma2_2_num + r[i,k]*t[g,h]*(y[i,g]-nu[k,h])^2
            sigma2_2_den = sigma2_2_den + r[i,k]*t[g,h]
          }
        }
      }
    }
    sigma2_2 = sigma2_2_num/sigma2_2_den
  
      
    # free energy
    Free_E = 0
    for(i in 1:N){
      for(k in 1:P){
        for(j in min((i+1),N):N){
          for(l in 1:P){
            Free_E = Free_E + r[i,k]*r[j,l]*(-((x[i,j] - mu[k,l])^2)/(2*sigma2_1) - log(sqrt(2*pi*sigma2_1)))
          }
        }
        Free_E = Free_E + log(rho[k]^r[i,k]) - log(r[i,k]^r[i,k])
      }
    }
    for(g in 1:M){
      for(h in 1:Q){
        for(i in 1:N){
          for(k in 1:P){
            Free_E = Free_E + r[i,k]*t[g,h]*(-((y[i,g] - nu[k,h])^2)/(2*sigma2_2) - log(sqrt(2*pi*sigma2_2)))
          }
        }
        Free_E = Free_E + log(tau[h]^t[g,h]) - log(t[g,h]^t[g,h])
      }
    }
    
    if(abs((Free_E-Free_E0)/Free_E)<eps){
      break
    }
      Free_E0 = Free_E
      print(Free_E)
  }
  
  return(list("r"=r,"t"=t,"rho"=rho,"tau"=tau,"mu"=mu,"nu"=nu,"sig1"=sigma2_1,"sig2"=sigma2_2))
}
