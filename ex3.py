#exercicio 3 - final 

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import matplotlib
matplotlib.rc('font', size=18)
matplotlib.rc('font', family='Arial')
file=open('n3.dat','w')
N = 501 
dt = 5.e-5 
L = float(2) 
nsteps = 1000 
dx = L/(N-1) 
nplot = 20 
r = dt/dx**2
A = np.zeros((N,N))
B = np.zeros((N,N))
for i in range(N):
    if i==0:
        A[i,:] = [1+r if j==0 else -r if j==1 else 0 for j in range(N)]
    elif i==N-1:
        A[i,:] = [1+r if j==N-1 else -r if j==N-2 else 0 for j in range(N)]
    else:
        A[i,:] = [-r if j==i-1 or j==i+1 else 1+2*r if j==i else 0 for j in range(N)]
x = np.linspace(0,1,N)
n=3
k=2*np.pi*10*n/(2.0*L)
u = np.asarray([np.cos(k*xx) for xx in x])
Amp = []
fig = plt.figure()
plt.plot(x,u,linewidth=2)
filename = 'foo000.jpg';
fig.set_tight_layout(True);
plt.xlabel("x")
plt.ylabel("u")
plt.title("t = 0")
plt.savefig(filename)
plt.clf()
c = 0
time=0.
tempo = []
for j in range(nsteps):
    u[:] = np.linalg.solve(A,u)
    Amp.append(max(u))
    time = time + dt
    tempo.append(time)
    print(tempo[j],Amp[j],file=file)
    print(j)
    bb = B.dot(u[:]) 
    if(j%nplot==0):
        plt.plot(x,u,linewidth=2)
        plt.ylim([-1,1])
        filename = 'foo' + str(c+1).zfill(3) + '.jpg';
        plt.xlabel("x")
        plt.ylabel("u")
        plt.title("t = %2.2f"%(dt*(j+1)))
        plt.savefig(filename)
        plt.clf()
        c += 1
plt.plot(tempo,Amp,linewidth=2)
filename = 'foo' + str(c+1).zfill(3) + '.jpg';
plt.xlabel("t")
plt.ylabel("A(t)")
plt.title("Amplitude em função do tempo")
plt.savefig("amp.jpg")
plt.clf()
os.system("ffmpeg -r 5 -y -i 'foo%03d.jpg' ex3-1.m4v")
os.system("rm -f *.jpg")
