# all the libraries required

library(e1071)
library(bikm1)
library(blockmodels)
library(fcd)
library(ggplot2)
library(lattice)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SBM_LBM_gen.R")
source("init_clust_spectral.R")
source("easy_init.R")
source("SBM_EMV.R")
source("LBM_EMV.R")
source("SLBM_EMV.R")
source("test.R")
source("ICL.R")
source("EMV_Select_model_SLBM.R")



# Generate and save data

N = 300
M = 300
P = 4
Q = 4
set.seed(123)

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
    save(list = list_name, file = paste(folder,file_name,sep="/"))
  }
}





# Calcul misc,ARI/CARI

calcul_SBM_erreur = function(x,y,z,w){
  
  # initialisation
  init = init_cs(x,y,P,Q)
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

# calcul_easySBM_erreur = function(x,y,z,w){
# 
#   model_BM = BM_gaussian("SBM_sym", x)
#   model_BM$estimate()
#   r = model_BM$memberships[[P]]$Z
#   
#   return(list("misc"=misc(r,z),"ARI"=get_ARI(r,z)))
# }

calcul_LBM_erreur = function(x,y,z,w){
  
  # initialisation
  init = init_cs(x,y,P,Q)
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
  init = init_cs(x,y,P,Q)
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




# save results for SBM

Tab_misc = matrix(0,50,5)
Tab_ARI = matrix(0,50,5)
for(i in 1:5){
  folder = paste("Data/test1/SLBM_theta",as.character(i))
  for(k in 1:50){
    file_name = paste("data",as.character(k),sep="")
    load(file = paste(folder,file_name,sep="/"))
    
    print(k)
    ERR = calcul_SBM_erreur(x,y,z,w)
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
save(list = list_name, file = paste(folder,file_name,sep="/"))

folder = paste("Result/test1/SBM0")
file_name = "result"
load(paste(folder,file_name,sep="/"))
data_plot = data.frame(sigma = rep(Tsigma,each=50), misc=matrix(Tab_misc,50*5,1), ARI=matrix(Tab_ARI,50*5,1))
ggplot(data_plot, aes(x = factor(sigma), y = misc)) + geom_boxplot()
ggplot(data_plot, aes(x = factor(sigma), y = ARI)) + geom_boxplot() 



# save results for LBM

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
save(list = list_name, file = paste(folder,file_name,sep="/"))


folder = paste("Result/test2_LBM_theta")
file_name = "result"
load(paste(folder,file_name,sep="/"))
data_plot = data.frame(sigma = rep(Tsigma,each=50), misc=matrix(Tab_misc,50*5,1), CARI=matrix(Tab_CARI,50*5,1))
ggplot(data_plot, aes(x = factor(sigma), y = misc)) + geom_boxplot()
ggplot(data_plot, aes(x = factor(sigma), y = CARI)) + geom_boxplot()



# save results for SLBM

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
save(list = list_name, file = paste(folder,file_name,sep="/"))


folder = paste("Result/test1/SLBM")
file_name = "result"
load(paste(folder,file_name,sep="/"))
data_plot = data.frame(sigma = rep(Tsigma,each=50), misc=matrix(Tab_misc,50*5,1), CARI=matrix(Tab_CARI,50*5,1))
ggplot(data_plot, aes(x = factor(sigma), y = misc)) + geom_boxplot()
ggplot(data_plot, aes(x = factor(sigma), y = CARI)) + geom_boxplot()

rep(Tsigma,each=50)







# load data from GLV and ssave results


N = M = 100
P = 2
Q = 3

file = "sigma04"
#file = ""

X = read.table(paste("Data/TestLV/Xfile",file,".txt",sep=""), header = FALSE, sep = "", dec = ".")
x = matrix(0,N,N)
for(i in 1:N){
  for(j in 1:N){
    x[i,j] = X[i,j]
  }
}
Y = read.table(paste("Data/TestLV/Yfile",file,".txt",sep=""), header = FALSE, sep = "", dec = ".")
y = matrix(0,N,M)
for(i in 1:N){
  for(g in 1:M){
    y[i,g] = Y[i,g]
  }
}

init = init_cs(x,y,P,Q)
#init = easy_init(x,y,z,w)
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

Z = read.table(paste("Data/TestLV/zlabel",file,".txt",sep=""))
W = read.table(paste("Data/TestLV/wlabel",file,".txt",sep=""))

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
table(apply(z,1,which.max),apply(r,1,which.max))

z_ind = apply(z,1,which.max)
w_ind = apply(w,1,which.max)
levelplot(t(x[order(z_ind),order(z_ind)]),col.regions=heat.colors(100))
levelplot(t(y[order(z_ind),order(w_ind)]),col.regions=heat.colors(100))

z2 = apply(r,1,which.max)
z3 = rep(0,N)
for(i in 1:N){
  if(z2[i] == 1){
    z3[i] = 2
  }
  if(z2[i] == 2){
    z3[i] = 1
  }
}
w2 = apply(t,1,which.max)
w3 = rep(0,M)
for(i in 1:M){
  if(w2[i] == 1){
    w3[i] = 3
  }
  if(w2[i] == 2){
    w3[i] = 2
  }
  if(w2[i] == 3){
    w3[i] = 1
  }
}

levelplot(t(y),col.regions=heat.colors(100))
levelplot(t(y[order(z3),order(w3)]),col.regions=heat.colors(100))

write.table(z,paste("Result/TestLV/zlabel",file,".txt",sep=""))
write.table(w,paste("Result/TestLV/wlabel",file,".txt",sep=""))





# EMV + Model selection

N = 100
M = 100

for(P in 4:8){
  folder = paste("Data/testICL/SLBM",as.character(P))
  if(! file.exists(folder)) {
    dir.create(folder)
  }
  for(k in 1:50){
    Q = P
    SLBM = SLBM_gen2(M,N,P,Q,epsilon=0.1,sigma2=0.1)
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
    save(list = list_name, file = paste(folder,file_name,sep="/"))
  }
}


ICL_SLBM_erreur = function(x,y){
  
  # EMV
  res_EMV = EMV_Select_model(x,y)
  r = res_EMV$r
  t = res_EMV$t
  rho = res_EMV$rho
  tau = res_EMV$tau
  mu = res_EMV$mu
  nu = res_EMV$nu
  sigma2_1 = res_EMV$sig1
  sigma2_2 = res_EMV$sig2
  
  return(list("P"=dim(r)[2],"Q"=dim(t)[2]))
}


Tab_P = matrix(0,4,50)
Tab_Q = matrix(0,4,50)
for(P in 4:8){
  folder = paste("Data/testICL/SLBM",as.character(P))
  for(k in 1:50){
    file_name = paste("data",as.character(k),sep="")
    load(file = paste(folder,file_name,sep="/"))
    
    ERR = ICL_SLBM_erreur(x,y)
    Tab_P[P-3,k] = ERR$P
    Tab_Q[P-3,k] = ERR$Q
  }
}
folder = paste("Result/testICL/SLBM")
if(! file.exists(folder)){
  dir.create(folder)
}
file_name = "result"
list_name = c("Tab_P","Tab_Q")
save(list = list_name, file = paste(folder,file_name,sep="/"))


Sol = EMV_Select_model(x,y)
r = Sol$r
t = Sol$t
dim(r)[2]
dim(t)[2]
