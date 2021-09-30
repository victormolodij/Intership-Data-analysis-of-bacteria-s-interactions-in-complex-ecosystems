import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import csv

N = 6
J = 3
K = 200
T = 20

mu = [7.5,2.6,2.5]*2
A = [[-2,-5,-0.5]*2,[-0.5,-1,-1.2]*2,[-1,-0.5,-1]*2]*2

X = np.linspace(0,T,K)


Y0 = [[1,1,1,0.5,0.5,0.5],[1,1,1,0.5,0.5,0.5],[1,1,1,0.5,0.5,0.5]]
Y = [[]]*J

def f(y,t):
  return((y * (mu + np.dot(A,y))))

for i in range(J):
    Y[i] = scipy.integrate.odeint(f,Y0[i],X)
    
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(X,Y[i])
    ax.legend(["x1","x2","x3","x4","x5","x6"])
    plt.show()
    


ax = plt.subplot(111)
ax.plot(X,Y[0][:,0])
ax.plot(X,Y[0][:,3])
ax.legend(["x1","x4"])
plt.show()


"""
Y = np.zeros((K*J,N))
Y[0] = Y0[0]
for i in range(1,K):
    Y[i] = Y[i-1] + (X[i] - X[i-1])*f(Y[i-1],X[i-1])

plt.figure()
plt.plot(X,Y)
plt.show()
"""

Z = np.zeros((N,K*J))
for n in range(N):
    for j in range(J):
        for k in range(K):
            Z[n,K*j + k] = Y[j][k][n]


csvfile = "data.csv"

with open(csvfile, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(Z)
