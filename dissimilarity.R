# table_especes est une table dont les lignes sont les espèces bactériennes
# les colonnes sont les temps de mesures ordonnés par individus (souris) et temps de mesure

#npt est un vecteur ordonné par individu qui pour chaque individu contient le nombre de temps de mesure

# pairwise cosine similarity between 2 species (species ix[1] and species ix[2])
# first computed in each individual mouse (0 if one in the pair is absent in the mouse) 
# weighted mean over the mice where both are present, weights are the number of observation in each mouse
cos.simTp = function(ix) 
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
    pN=sqrt(sum(pA^2)*sum(pB^2))
    if (pN>0) {
      S = S + npt[i]*sum(pA*pB)/pN
      T = T + npt[i]
    }
    offset=b+1
  }
  if (T==0) {
    T=1
    #print("T==0")
    } #this is not good T=0 just means that these species are never found in the same mouse. What should we do?
  return(S/T)
}