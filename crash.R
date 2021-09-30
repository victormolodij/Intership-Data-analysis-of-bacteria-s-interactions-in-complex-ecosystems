SLBM = SLBM_gen2(M,N,P,Q,epsilon=0.1,sigma2=5)
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



image(x)
z_ind = apply(z,1,which.max)
image(x[order(z_ind),order(z_ind)])
#image(y)
#w_ind = apply(w,1,which.max)
#image(y[order(z_ind),order(w_ind)])

model_BM = BM_gaussian("SBM_sym", x)
model_BM$estimate()
r = model_BM$memberships[[P]]$Z
misc(r,z)

# initialisation
#init = easy_init(x,y,z,w)
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
#res_EMV = SLBM_EMV(x,y,r,t,rho,tau,mu,nu,sigma2_1,sigma2_2)
res_EMV = SBM_EMV(x,r,rho,mu,sigma2_1)
#res_EMV = LBM_EMV(y,r,t,rho,tau,nu,sigma2_2)
r = res_EMV$r
t = res_EMV$t


misc(r,z)
misc(t,w)
