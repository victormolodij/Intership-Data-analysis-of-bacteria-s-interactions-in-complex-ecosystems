ICL = function(N,M,P,Q){
  return(0)
}

a = 1
b = 2
a=b=a



selec_model(z,w){
  
  ICL0 = -Inf
  
  while(ICL_z_plus > ICL0 || ICL_w_plus > ICL0){
    
    z_index = sample(1:P, 1)
    z_plus = cbind(z,matrix(0,N,1))
    for(i in 1:N){
      rd = sample(0:1,1)
      if(rd){
        z_plus[i,P+1] = z_plus[i,z_index]
        z_plus[i,z_index] = 0
      }
    }
    
    w_index = sample(1:Q, 1)
    w_plus = cbind(t,matrix(0,M,1))
    for(g in 1:M){
      rd = sample(0:1,1)
      if(rd){
        w_plus[g,Q+1] = w_plus[g,w_index]
        w_plus[g,w_index] = 0
      }
    }


    # ICL
    ICL_r_plus = ICL()
    ICL_t_plus = ICL()
    if(ICL_r_plus >= ICL_t_plus){
      P = P+1
      ICL0 = ICL_r_plus
    }
    else{
      Q = Q+1
      ICL0 = ICL_t_plus
    }
  }
}


