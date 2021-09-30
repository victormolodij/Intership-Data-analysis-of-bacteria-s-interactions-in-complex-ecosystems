init_cs = function(x,y){
  CS_x = spectral.clustering(x, normalised = TRUE, K = P)
  CS_y = spectral.clustering(y, normalised = TRUE, K = Q)
  mem_x = CS_x[1:N]
  mem_y = CS_y[1:N]
  
  r = t(matrix(as.numeric(rep(mem_x,each=P) == 1:P),P,N))
  t = t(matrix(as.numeric(rep(mem_y,each=Q) == 1:Q),Q,M))
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

  return(list("r"=r,"t"=t,"rho"=rho,"tau"=tau,"mu"=mu,"nu"=nu,"sig1"=sigma2_1,"sig2"=sigma2_2))
}
