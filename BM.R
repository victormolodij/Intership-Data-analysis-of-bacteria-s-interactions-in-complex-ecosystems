model_BM = BM_gaussian("SBM_sym", CTp)
model_BM$estimate()
Q = which.max(model_BM$ICL)
Z = model_BM$memberships[[Q]]$Z
Z = apply(Z,1,which.max)
image(CTp)
image(CTp[order(Z),order(Z)])









model_BM = BM_gaussian("SBM_sym", x, explore_min = P, explore_max = P, exploration_factor = 1)
model_BM$estimate()
r = model_BM$memberships[[P]]$Z
misc(r,z)
