# Model selection and EMV on SLBM

EMV_Select_model = function(x,y,P_min=2,Q_min=2){
  
  N = dim(x)[1]
  M = dim(y)[2]
  P = P_min
  Q = Q_min
  
  init = init_cs(x,y,P,Q)
  res_EMV = SLBM_EMV(x,y,init$r,init$t,init$rho,init$tau,init$mu,init$nu,init$sig1,init$sig2)
  r = res_EMV$r
  t = res_EMV$t
  ICL = ICL_asp(x,y,r,t,init$rho,init$tau,init$mu,init$nu,init$sig1,init$sig2)
  
  while(TRUE){
    ICL_P0 = -Inf
    for(k in 1:P){
      r_plus = cbind(r,matrix(0,N,1))
      smpl = sample(1:N,N/2)
      r_plus[smpl,k] = 0
      r_plus[smpl,P+1] = r[smpl,k]
      
      init = easy_init(x,y,r_plus,t)
      res_EMV = SLBM_EMV(x,y,init$r,init$t,init$rho,init$tau,init$mu,init$nu,init$sig1,init$sig2)
      r_plus = res_EMV$r
      ICL_P = ICL_asp(x,y,r_plus,res_EMV$t,res_EMV$rho,res_EMV$tau,res_EMV$mu,res_EMV$nu,res_EMV$sig1,res_EMV$sig2)
      if(ICL_P>ICL_P0){
        ICL_P0 = ICL_P
        r_save = r_plus
      }
    }
    
    ICL_Q0 = -Inf
    for(h in 1:Q){
      t_plus = cbind(t,matrix(0,M,1))
      smpl = sample(1:M,M/2)
      t_plus[smpl,h] = 0
      t_plus[smpl,Q+1] = t[smpl,h]
      
      init = easy_init(x,y,r,t_plus)
      res_EMV = SLBM_EMV(x,y,init$r,init$t,init$rho,init$tau,init$mu,init$nu,init$sig1,init$sig2)
      t_plus = res_EMV$t
      ICL_Q = ICL_asp(x,y,res_EMV$r,t_plus,res_EMV$rho,res_EMV$tau,res_EMV$mu,res_EMV$nu,res_EMV$sig1,res_EMV$sig2)
      if(ICL_Q>ICL_Q0){
        ICL_Q0 = ICL_Q
        t_save = t_plus
      }
    }
    if(ICL_P0>ICL | ICL_Q0>ICL){
      if(ICL_P0>ICL_Q0){
        r = r_save
        ICL = ICL_P0
        P = P+1
      }
      else{
        t = t_save
        ICL = ICL_Q0
        Q = Q+1
      }
      print(c("P=",P))
      print(c("Q=",Q))
    }
    else{
      break()
    }
  }
  return(list("r"=r,"t"=t,"rho"=rho,"tau"=tau,"mu"=mu,"nu"=nu,"sig1"=sigma2_1,"sig2"=sigma2_2))
}
