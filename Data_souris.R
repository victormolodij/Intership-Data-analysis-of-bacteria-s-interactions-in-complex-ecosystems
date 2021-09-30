Data = read.csv("Data/Data_souris.csv")

npt = table(Data$Individu)
npt = c(npt)

table_especes = Data[-1:-2]
table_especes = t(table_especes)

#table_especes = table_especes[1:50,]
dim(table_especes)



# Similarity

n <- nrow(table_especes) 
cmb <- expand.grid(i=1:n, j=1:n) 

CTp <- matrix(apply(cmb,1,cos.simTp),n,n)
CosinusTp <- as.dist(m=(1-CTp))