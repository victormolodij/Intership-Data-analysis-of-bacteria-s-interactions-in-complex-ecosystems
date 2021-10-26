# Calcul of ICL asymptotique

ICL_asp = function(x,y,r,t,rho,tau,mu,nu,sigma2_1,simga2_2){
  
  P = dim(r)[2]
  N = dim(r)[1]
  Q = dim(t)[2]
  M = dim(t)[1]
  
  ICL_value = 0
  z = apply(r,1,which.max)
  w = apply(t,1,which.max)
  for(i in 1:N){
    for(j in 1:N){
      ICL_value = ICL_value + (-((x[i,j] - mu[z[i],z[j]])^2)/(2*sigma2_1) - log(sqrt(2*pi*sigma2_1)))
    }
    ICL_value = ICL_value + rho[z[i]]
  }
  for(g in 1:M){
    for(i in 1:N){
      ICL_value = ICL_value + (-((y[i,g] - nu[z[i],w[g]])^2)/(2*sigma2_2) - log(sqrt(2*pi*sigma2_2)))
    }
  }
  ICL_value = ICL_value + tau[w[g]]
  
  ICL_value = ICL_value - (P*(P+1)/4)*log(N*(N-1)/2) - (P*Q/2)*log(N*M) - (P/2)*log(N) - (Q/2)*log(M)
  return(ICL_value)
}
