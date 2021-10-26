# Misclassified rate

misc = function(r,z){
  N = dim(z)[1]
  P = dim(z)[2]
  Permu = permutations(P)
  Tmem = apply(r,1,which.max)
  Tmem0 = matrix(0,N,P)
  for(i in 1:N){
    Tmem0[i,Tmem[i]] = 1
  }
  
  TP = rep(0,dim(Permu)[1])
  for(i in 1:N){
    for(k in 1:P){
      for(sig in 1:dim(Permu)[1]){
        if(z[i,k]==0 && Tmem0[i,Permu[sig,k]]==1){
          TP[sig] = TP[sig] + 1
        }
      }
    }
  }
  return(min(TP)/N)
}


# ARI/CARI

get_ARI = function(r,z){
  N = dim(z)[1]
  P = dim(z)[2]
  Tmem = apply(r,1,which.max)
  Tmem0 = matrix(0,N,P)
  for(i in 1:N){
    Tmem0[i,Tmem[i]] = 1
  }
  return(ARI(apply(z,1,which.max), apply(Tmem0,1,which.max)))
}

get_CARI = function(r,t,z,w){
  N = dim(z)[1]
  P = dim(z)[2]
  Tmem1 = apply(r,1,which.max)
  Tmem1_0 = matrix(0,N,P)
  for(i in 1:N){
    Tmem1_0[i,Tmem1[i]] = 1
  }
  
  M = dim(w)[1]
  Q = dim(w)[2]
  Tmem2 = apply(t,1,which.max)
  Tmem2_0 = matrix(0,M,Q)
  for(i in 1:M){
    Tmem2_0[i,Tmem2[i]] = 1
  }
  return(CARI(apply(z,1,which.max),apply(w,1,which.max),apply(Tmem1_0,1,which.max),apply(Tmem2_0,1,which.max)))
}
