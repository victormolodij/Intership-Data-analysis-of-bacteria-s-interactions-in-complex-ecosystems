# Load data generated with GLV.py

Y = read.csv("data.csv", header= FALSE)

# number of groupe in GLV
J = 3


table_especes = NULL
for (i in 1:(10*J)){
  table_especes = cbind(table_especes,Y[,20*i])
}
npt = rep(10,J)


b = rnorm(N*K*J, mean = 0, sd = 1)


matplot(t(table_especes[1:3,1:10]))