# generate synthetic data with SLBM model


SLBM_gen = function(M,N,P,Q){
  epsilon1 = 0.4
  mu = (1-epsilon1)*matrix(1,P,P)
  diag(mu) = epsilon1
  epsilon2 = 0.4
  nu = (1-epsilon2)*matrix(1,P,Q)
  diag(nu) = epsilon2
  
  sigma2_1 = 0.1
  sigma2_2 = 0.1
  rho = rep(1/P,P)
  tau = rep(1/Q,Q)
  z = t(rmultinom(N, 1, rho))
  w = t(rmultinom(M, 1, tau))
  x = matrix(rnorm(N*N,0,sqrt(sigma2_1)),N,N) + z%*%mu%*%t(z)
  x[lower.tri(x)] = t(x)[lower.tri(x)]
  diag(x) = 0
  y = matrix(rnorm(N*M,0,sqrt(sigma2_2)),N,M) + z%*%nu%*%t(w)
  
  return(list("x"=x,"y"=y,"z"=z,"w"=w,"rho"=rho,"tau"=tau,"mu"=mu,"nu"=nu,"sig1"=sigma2_1,"sig2"=sigma2_2))
}



SLBM_gen1 = function(M,N,P,Q,epsilon,sigma2){
  epsilon1 = epsilon
  mu = epsilon1*matrix(1,P,P)
  diag(mu) = 1-epsilon1
  epsilon2 = epsilon
  nu = epsilon2*matrix(1,P,Q)
  diag(nu) = 1-epsilon2
  sigma2_1 = sigma2
  sigma2_2 = sigma2
  rho = rep(1/P,P)
  tau = rep(1/Q,Q)
  z = t(rmultinom(N, 1, rho))
  w = t(rmultinom(M, 1, tau))
  x = matrix(rnorm(N*N,0,sqrt(sigma2_1)),N,N) + z%*%mu%*%t(z)
  x[lower.tri(x)] = t(x)[lower.tri(x)]
  diag(x) = 0
  y = matrix(rnorm(N*M,0,sqrt(sigma2_2)),N,M) + z%*%nu%*%t(w)
  
  return(list("x"=x,"y"=y,"z"=z,"w"=w,"rho"=rho,"tau"=tau,"mu"=mu,"nu"=nu,"sig1"=sigma2_1,"sig2"=sigma2_2))
}



SLBM_gen2 = function(M,N,P,Q,epsilon,sigma2){
  epsilon1 = epsilon
  mu = matrix(0,P,P)
  for(k in 1:P){
    D = rbind(cbind(matrix(0,k,P-k),diag(P-k,k,k)),matrix(0,P-k,P))
    mu = mu + epsilon1*D
  }
  mu[lower.tri(mu)] = t(mu)[lower.tri(mu)]
  diag(mu) = 1-epsilon1
  epsilon2 = epsilon
  nu = t(epsilon2*matrix(seq(1,P),Q,P))
  diag(nu) = 1-epsilon2
  sigma2_1 = sigma2
  sigma2_2 = sigma2
  rho = rep(1/P,P)
  tau = rep(1/Q,Q)
  z = t(rmultinom(N, 1, rho))
  w = t(rmultinom(M, 1, tau))
  x = matrix(rnorm(N*N,0,sqrt(sigma2_1)),N,N) + z%*%mu%*%t(z)
  x[lower.tri(x)] = t(x)[lower.tri(x)]
  diag(x) = 0
  y = matrix(rnorm(N*M,0,sqrt(sigma2_2)),N,M) + z%*%nu%*%t(w)
  
  return(list("x"=x,"y"=y,"z"=z,"w"=w,"rho"=rho,"tau"=tau,"mu"=mu,"nu"=nu,"sig1"=sigma2_1,"sig2"=sigma2_2))
}





# test

# N = 200
# M = 200
# P = 4
# Q = 4
# SLBM = SLBM_gen(M,N,P,Q)
# x = SLBM$x
# y = SLBM$y
# z = SLBM$z
# w = SLBM$w
# rho = SLBM$rho
# tau = SLBM$tau
# mu = SLBM$mu
# nu = SLBM$nu
# sigma2_1 = SLBM$sig1
# sigma2_2 = SLBM$sig2
# 
# k = 2
# l = 2
# Gk = (z[,k]*(1:N))[z[,k] != 0]
# Gl = (z[,l]*(1:N))[z[,l] != 0]
# mean((x[Gk,Gl])[lower.tri(x[Gk,Gl])])
# mu[k,l]
# 
# k = 2
# h = 2
# Gg = (z[,k]*(1:N))[z[,k] != 0]
# Gh = (w[,h]*(1:M))[w[,h] != 0]
# mean(y[Gg,Gh])
# nu[k,h]