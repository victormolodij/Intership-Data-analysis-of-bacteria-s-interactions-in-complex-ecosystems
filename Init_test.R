# test



r0 = z
r = z
r[1:50,1] = r0[1:50,2]
r[1:50,2] = r0[1:50,1]
rho = (1/P)*matrix(1,P,1)
mu = mu0 + 0.1

t0 = w
t = w
t[1:25,1] = t0[1:25,2]
t[1:25,2] = t0[1:25,1]
tau = (1/Q)*matrix(1,Q,1)
nu = nu0 + 0.1