library(e1071)
library(bikm1)
library(blockmodels)
library(fcd)
library(ggplot2)



N = 300
M = 300
P = 4
Q = 4
set.seed(123)



# Generate and save data

Tsigma = c(0.1,1,5,10,20)
for(i in 1:5){
  folder = paste("Data/test1/SLBM_theta",as.character(i))
  if(! file.exists(folder)) {
    dir.create(folder)
  }
  for(k in 1:50){
    SLBM = SLBM_gen2(M,N,P,Q,epsilon=0.1,sigma2=Tsigma[i])
    x = SLBM$x
    y = SLBM$y
    z = SLBM$z
    w = SLBM$w
    rho = SLBM$rho
    tau = SLBM$tau
    mu = SLBM$mu
    nu = SLBM$nu
    sigma2_1 = SLBM$sig1
    sigma2_2 = SLBM$sig2
    
    file_name = paste("data",as.character(k),sep="")
    list_name = c("x","y","z","w","rho","tau","mu","nu","sigma2_1","sigma2_2")
    #save(list = list_name, file = paste(folder,file_name,sep="/"))
  }
}





# Calcul misc,ARI/CARI

calcul_SBM_erreur = function(x,y,z,w){
  
  # initialisation
  init = init_cs(x,y)
  #init = easy_init(x,y,z,w)
  r = init$r
  t = init$t
  rho = init$rho
  tau = init$tau
  mu = init$mu
  nu = init$nu
  sigma2_1 = init$sig1
  sigma2_2 = init$sig2
  
  # EMV
  res_EMV = SBM_EMV(x,r,rho,mu,sigma2_1)
  r = res_EMV$r
  rho = res_EMV$rho
  mu = res_EMV$mu
  sigma2_1 = res_EMV$sig1
  
  return(list("misc"=misc(r,z),"ARI"=get_ARI(r,z)))
}

calcul_easySBM_erreur = function(x,y,z,w){

  model_BM = BM_gaussian("SBM_sym", x)
  model_BM$estimate()
  r = model_BM$memberships[[P]]$Z
  
  return(list("misc"=misc(r,z),"ARI"=get_ARI(r,z)))
}

calcul_LBM_erreur = function(x,y,z,w){
  
  # initialisation
  init = init_cs(x,y)
  r = init$r
  t = init$t
  rho = init$rho
  tau = init$tau
  mu = init$mu
  nu = init$nu
  sigma2_1 = init$sig1
  sigma2_2 = init$sig2
  
  # EMV
  res_EMV = LBM_EMV(y,r,t,rho,tau,nu,sigma2_2)
  r = res_EMV$r
  t = res_EMV$t
  rho = res_EMV$rho
  tau = res_EMV$tau
  nu = res_EMV$nu
  sigma2_2 = res_EMV$sig2
  
  return(list("misc"= 0.5*(misc(r,z)+misc(t,w)) ,"CARI"=get_CARI(r,t,z,w)))
}

calcul_SLBM_erreur = function(x,y,z,w){
  
  # initialisation
  #init = init_cs(x,y)
  init = easy_init(x,y,z,w)
  r = init$r
  t = init$t
  rho = init$rho
  tau = init$tau
  mu = init$mu
  nu = init$nu
  sigma2_1 = init$sig1
  sigma2_2 = init$sig2
  
  # EMV
  res_EMV = SLBM_EMV(x,y,r,t,rho,tau,mu,nu,sigma2_1,sigma2_2)
  r = res_EMV$r
  t = res_EMV$t
  rho = res_EMV$rho
  tau = res_EMV$tau
  mu = res_EMV$mu
  nu = res_EMV$nu
  sigma2_1 = res_EMV$sig1
  sigma2_2 = res_EMV$sig2
  
  return(list("misc"= 0.5*(misc(r,z)+misc(t,w)) ,"CARI"=get_CARI(r,t,z,w)))
}



Tab_misc = matrix(0,50,5)
Tab_ARI = matrix(0,50,5)
for(i in 1:5){
  folder = paste("Data/test1/SLBM_theta",as.character(i))
  for(k in 1:50){
    file_name = paste("data",as.character(k),sep="")
    load(file = paste(folder,file_name,sep="/"))
    
    ERR = calcul_easySBM_erreur(x,y,z,w)
    Tab_misc[k,i] = ERR$misc
    Tab_ARI[k,i] = ERR$ARI$ari
  }
}


folder = paste("Result/test1/SBM")
if(! file.exists(folder)) {
  dir.create(folder)
}
file_name = "result"
list_name = c("Tab_misc","Tab_ARI")
#save(list = list_name, file = paste(folder,file_name,sep="/"))

folder = paste("Result/test1/SBM0")
file_name = "result"
load(paste(folder,file_name,sep="/"))
data_plot = data.frame(sigma = rep(Tsigma,each=50), misc=matrix(Tab_misc,50*5,1), ARI=matrix(Tab_ARI,50*5,1))
ggplot(data_plot, aes(x = factor(sigma), y = misc)) + geom_boxplot()
ggplot(data_plot, aes(x = factor(sigma), y = ARI)) + geom_boxplot() 





Tab_misc = matrix(0,50,5)
Tab_CARI = matrix(0,50,5)
for(i in 1:5){
  folder = paste("Data/test_SLBM_theta",as.character(i))
  for(k in 1:50){
    file_name = paste("data",as.character(k),sep="")
    load(file = paste(folder,file_name,sep="/"))
    
    ERR = calcul_LBM_erreur(x,y,z,w)
    Tab_misc[k,i] = ERR$misc
    Tab_CARI[k,i] = ERR$CARI$cari
  }
}

folder = paste("Result/test2/LBM_theta")
if(! file.exists(folder)) {
  dir.create(folder)
}
file_name = "result"
list_name = c("Tab_misc","Tab_CARI")
#save(list = list_name, file = paste(folder,file_name,sep="/"))


folder = paste("Result/test2_LBM_theta")
file_name = "result"
load(paste(folder,file_name,sep="/"))
data_plot = data.frame(sigma = rep(Tsigma,each=50), misc=matrix(Tab_misc,50*5,1), CARI=matrix(Tab_CARI,50*5,1))
ggplot(data_plot, aes(x = factor(sigma), y = misc)) + geom_boxplot()
ggplot(data_plot, aes(x = factor(sigma), y = CARI)) + geom_boxplot()





Tab_misc = matrix(0,50,5)
Tab_CARI = matrix(0,50,5)
for(i in 1:5){
  folder = paste("Data/test1/SLBM_theta",as.character(i))
  for(k in 1:50){
    file_name = paste("data",as.character(k),sep="")
    load(file = paste(folder,file_name,sep="/"))
    
    ERR = calcul_SLBM_erreur(x,y,z,w)
    Tab_misc[k,i] = ERR$misc
    Tab_CARI[k,i] = ERR$CARI$cari
  }
}

folder = paste("Result/test1/SLBM")
if(! file.exists(folder)) {
  dir.create(folder)
}
file_name = "result"
list_name = c("Tab_misc","Tab_CARI")
#save(list = list_name, file = paste(folder,file_name,sep="/"))


folder = paste("Result/test1/SLBM")
file_name = "result"
load(paste(folder,file_name,sep="/"))
data_plot = data.frame(sigma = rep(Tsigma,each=50), misc=matrix(Tab_misc,50*5,1), CARI=matrix(Tab_CARI,50*5,1))
ggplot(data_plot, aes(x = factor(sigma), y = misc)) + geom_boxplot()
ggplot(data_plot, aes(x = factor(sigma), y = CARI)) + geom_boxplot()

rep(Tsigma,each=50)





N = M = 100
P = 2
Q = 3

X = read.table("Data/TestLV/Xfile.txt")
x = matrix(0,N,N)
for(i in 1:N){
  for(j in 1:N){
    x[i,j] = X[i,j]
  }
}
Y = read.table("Data/TestLV/Yfile.txt", header = FALSE, sep = "", dec = ".")
y = matrix(0,N,M)
for(i in 1:N){
  for(g in 1:M){
    y[i,g] = Y[i,g]
  }
}

#init = init_cs(x,y)
init = easy_init(x,y,z,w)
r = init$r
t = init$t
rho = init$rho
tau = init$tau
mu = init$mu
nu = init$nu
sigma2_1 = init$sig1
sigma2_2 = init$sig2

res_EMV = SLBM_EMV(x,y,r,t,rho,tau,mu,nu,sigma2_1,sigma2_2)
r = res_EMV$r
t = res_EMV$t
rho = res_EMV$rho
tau = res_EMV$tau
mu = res_EMV$mu
nu = res_EMV$nu
sigma2_1 = res_EMV$sig1
sigma2_2 = res_EMV$sig2

Z = read.table("Data/TestLV/zlabel.txt")
W = read.table("Data/TestLV/wlabel.txt")

#apply(t,1,which.max)

z = matrix(0,N,P)
for(i in 1:N){
  z[i,Z[i,1]] = 1
}
w = matrix(0,M,Q)
for(g in 1:M){
  w[g,W[g,1]] = 1
}

misc(r,z)
misc(t,w)

table(apply(w,1,which.max),apply(t,1,which.max))

z = apply(r,1,which.max)
w = apply(t,1,which.max)

write.table(z,"Result/TestLV/Zfile.txt")
write.table(w,"Result/TestLV/Wfile.txt")
