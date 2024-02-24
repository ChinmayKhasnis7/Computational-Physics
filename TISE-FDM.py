## Chinmay Khasanis ##

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from scipy import sparse
from scipy.sparse import linalg as slinalg

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

def shrod_matrix(x_0, x_N, N, numeg):
    x = np.linspace(x_0, x_N, N)
    h = x[1]-x[0]
    V = potential(x)
    
    # Tridiagonal Hamiltonian matrix
    H = sparse.eye(N, N, format = "lil") * 2
    for i in range(N - 1):
        H[i, i + 1] = -1
        H[i + 1, i] = -1

    k = 1/(2*h**2)  # Natural units and mass =1
    H = H *k
    for i in range(N):
        H[i, i] = H[i, i] + V[i]
    H = H.tocsc()  # compressed sparse column format format
    # Eigenvalues and Eigenvectors of the Hamiltonian matrix using SciPy solver
    [E_opt, psi] = slinalg.eigs(H, k = numeg, which = "SM")
    E_opt = np.real(E_opt)
    # normalize
    for i in range(numeg):
        psi[:, i] = normalize_psi(psi[:, i], x)
    return E_opt, psi, x

x_0 = -5; x_N = 5; N = 1000; numeg = 10
E_opt, psi_arr, x_axis =  shrod_matrix(x_0 = -5, x_N = 5, N = 1000, numeg = 10)

#### Plot ####

plt.figure(figsize=(8,10))
for i in range(9):
    density = np.conj(psi_arr[:,i])*psi_arr[:,i]
    plt.plot(x_axis, density + E_opt[i], label='$E_{%i} = %.4f$' %(i,E_opt[i]))
x_axis_v  = np.linspace(-3.5,3.5, 2*N)   
plt.plot(x_axis_v,potential(x_axis_v),label='$Potential$')
plt.yticks(fontsize=12); plt.xticks(fontsize=12); plt.legend(fontsize=10)
plt.xlabel('$x$',fontsize=15)
plt.ylabel('$|\psi(x) |^2$ (shifted up by $E$ for visualisation)',fontsize=14)
plt.title('Harmonic oscillator: Matrix method',fontsize=12)
plt.show()