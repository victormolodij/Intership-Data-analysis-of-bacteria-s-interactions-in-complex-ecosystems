# table_especes est une table dont les lignes sont les espèces bactériennes
# les colonnes sont les temps de mesures ordonnés par individus (souris) et temps de mesure

#npt est un vecteur ordonné par individu qui pour chaque individu contient le nombre de temps de mesure

# pairwise cosine similarity between 2 species (species ix[1] and species ix[2])
# first computed in each individual mouse (0 if one in the pair is absent in the mouse) 
# weighted mean over the mice where both are present, weights are the number of observation in each mouse
support = function(L1,L2,n) 
{
  R = c()
  for (i in 1:n){
    if (L1[i] == 0 || L2[i] == 0){
      R = c(R,FALSE)
    }
    else{
      R = c(R,TRUE)
    }
  }
  return(R)
}

log.simTp = function(ix) 
{
  A = table_especes[ix[1],]
  B = table_especes[ix[2],]
  S = 0
  T = 0
  offset=1
  for (i in 1:length(npt)) {
    b=offset+npt[i]-1
    pA=A[offset:b]
    pB=B[offset:b]
    Supp = support(pA,pB,npt[i])
    if(sum(Supp) == 0){
      #print("next")
      next
    }
    pA = pA[Supp]
    pB = pB[Supp]
    pA = log(pA)
    pB = log(pB)
    pN = (pA - mean(pA)) - (pB - mean(pB))
    pN = norm(pN,type="2")
    S = S + pN
    T = T + 1
    offset=b+1
  }
  if (T==0) {
    #print("T==0")
    return(-1)
    } #this is not good T=0 just means that these species are never found in the same mouse. What should we do?
  return(S/T)
}

###############################################################################
#build similarity matrix

n <- nrow(table_especes) 
cmb <- expand.grid(i=1:n, j=1:n) 

CTp_log <- matrix(apply(cmb,1,log.simTp),n,n)

sum(is.na(CTp_log))

sum(CTp_log != t(CTp_log))

A = table_especes[1,1:npt[1]]
B = table_especes[2,1:npt[1]]
sum(support(A,B,npt[1]))
support(A,B,npt[1])
