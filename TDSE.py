## Chinmay Khasanis ##

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy import integrate


#### Infinite box potential
def box_potential(x):
    L = 5
    v0 = 0 ; v_inf = 100
    v = np.zeros_like(x)
    i=0
    for pt in x:
        if (0 < pt < L):
            v[i] = v0
        else:
            v[i] = v_inf
        i=i+1
    return v

# psi*psi
def prob_density(psi):
    return np.conj(psi)*psi

#define gaussian packet
def gauss_packet(x,x0,s):
    a = 1.0 / (s*np.sqrt(np.pi)) 
    return np.sqrt(a)*np.exp(-(x-x0)**2/(2*s**2))*np.exp(1j*0.1*x)

## Plot initial wavefunction at t=0
L=  5 # Note Change the L in potential function also if you change here
h = 0.01; s = 0.1 ; x0 = 2.5 
x = np.arange(0, 5+h, h) 
    
# Packet at t = 0
Psi_0 = gauss_packet(x,x0,s)
V = box_potential(x)
plt.figure(figsize=(10,5))
plt.plot(x,V,label='$V(x)$',color = "r")
plt.plot(x, prob_density(Psi_0), "g", label=r"$|\psi(t,x)|^2$")
plt.yticks(fontsize=12); plt.xticks(fontsize=12)
plt.xlabel('$x$',fontsize=15); plt.ylabel('$V(x)$',fontsize=14); plt.ylim(0,10)
plt.title('Gaussian wave packet at $t=0$',fontsize=12)
plt.legend(fontsize=12);plt.grid();plt.show()


# use FDM to get Laplace operator
D2 = sparse.diags([1, -2, 1], [-1, 0, 1], shape=(x.size, x.size))/h**2
D2.toarray()*h**2
# time evoln function using Laplace operator for solving IVP
def psi_t(t, psi):
    return -1j * (-0.5*D2.dot(psi) + V*psi)

# initial and final time
t0 = 0.0 ; tf = 0.1   
 # time steps
dt = 0.005 
t_eval = np.arange(t0, tf, dt)

# Initial Value Problem using built in Runge-Kutta
sol = integrate.solve_ivp(psi_t, t_span = [t0, tf], y0 = Psi_0, t_eval = t_eval, method="RK23")


## Plot time evolution of wavefunction
plt.figure(figsize=(10,10))
for i, t in enumerate(sol.t):
    plt.plot(x, np.abs(sol.y[:,i])**2, label='t=%.2f'%t)     
plt.plot(x,V,label='$V(x)$',color = "r")
plt.yticks(fontsize=12); plt.xticks(fontsize=12)
plt.xlabel('$x$',fontsize=15); plt.ylabel('$V(x)$',fontsize=14); plt.ylim(0,7)
plt.title('Gaussian wave packet evolution (Before hitting)',fontsize=12)
plt.legend(fontsize=12);plt.grid();plt.show()