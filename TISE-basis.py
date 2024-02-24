## Chinmay Khasanis ##

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from scipy import linalg

# Potential (V(x))
def potential(x):
    v = x**2
    return v

# psi*psi
def prob_density(psi):
    return np.conj(psi)*psi

# Normalization <psi|psi> = 1 
def normalize_psi(psi, x):
   area = np.sqrt(scipy.integrate.simps(prob_density(psi), x))
   return psi/area

###### 1d SHO using particle in box basis functions
def H_matrix(L,N):
    H= np.zeros([N,N],dtype=np.complex_)
    pi2_4 = 4#*np.pi**2
    L2 = L**2
    for m in range(N):
        for n in range(N): 
            if m==n: 
                H[m,n] = n**2/L2 + L2/pi2_4
            else:
                H[m,n] = (L2*(-1)**(m+n))/(pi2_4*(m-n)**2)
    return H

def correct_eigns(E_opt_, C_):
    E_opt_ = np.real(E_opt_)
    indices = np.arange(1,len(E_opt_)+1,1)
    sorted_indices = [x for _, x in sorted(zip(E_opt_, indices))]
    E_opt= np.zeros_like(E_opt_) ; C= np.zeros_like(C_)
    i=0
    for ind in sorted_indices:
        E_opt[i]= E_opt_[ind-1]
        C[:,i]= C_[:,ind-1]
        i=i+1
    s = [i for i in range(len(E_opt)) if E_opt[i] > 0][0]
    E_opt = E_opt[s:]
    C = C[:,:len(E_opt_)-s]

    return E_opt, C

def fourier_basis(x,n,l):
    phi_n = np.exp(n*0.5*1j*x)/np.sqrt(l)
    return phi_n

def back_to_cartesian(x,L,N,E_opt,C):
    psi = np.zeros([len(x),len(E_opt)], dtype=np.complex_)
    for i in range(len(E_opt)):
        summ = 0 + 0j
        coefs = C[:,i]
        for n in range(N):
            phi_n = fourier_basis(x,n,L)
            summ = summ + coefs[n]*phi_n
        psi[:,i] = summ #normalize_psi(summ,x)
    return psi


L = 5 ; N=20 ; h = 0.01
x = np.arange(-L,L,h)
H = H_matrix(L,N)
[E_opt_, C_] = linalg.eig(H, left=False, right=True)
E_opt, C = correct_eigns(E_opt_, C_)
psi_arr = back_to_cartesian(x,L,N,E_opt,C)


#### Plot ####

x_axis = np.linspace(-L,L,psi_arr.shape[0])
plt.figure(figsize=(8,10))
for i in range(9):
    density = np.conj(psi_arr[:,i])*psi_arr[:,i]
    plt.plot(x_axis, density + E_opt[i], label='$E_{%i} = %.4f$' %(i,E_opt[i]))
x_axis_v  = np.linspace(-3.5,3.5, 2*N)   
plt.plot(x_axis_v,potential(x_axis_v),label='$Potential$')
plt.yticks(fontsize=12); plt.xticks(fontsize=12); plt.legend(fontsize=9)
plt.xlabel('$x$',fontsize=15)
plt.ylabel('$|\psi(x) |^2$ (shifted up by $E$ for visualisation)',fontsize=14)
plt.title('Harmonic oscillator: Basis Method',fontsize=12)
plt.show()