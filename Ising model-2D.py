## Chinmay Khasanis ##

import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand

###########################  Subroutines ##################################
def montecarlo(spin_matrix, Boltzmann, N_spins):
    for _ in range(N_spins*N_spins):

        # choose a point in lattice randomly
        i = np.random.randint(0, N_spins)
        j = np.random.randint(0, N_spins)

        # compute the change in energy due to sign flip 
        dE = energy_change(spin_matrix, N_spins, i, j)

        # accept or reject this choice of flip 
        s_ij = metropolis_choice(spin_matrix[i,j], dE, Boltzmann)
        spin_matrix[i,j] = s_ij
    return spin_matrix

# Metropolis sampling condition
def metropolis_choice(s_ij, dE, Boltzmann):
    if dE < 0:
        s_ij = -1*s_ij
    elif rand() < np.exp(-dE*Boltzmann):
        s_ij = -1*s_ij
    return s_ij

# compute the change in energy due to sign flip 
def energy_change(spin_matrix, N_spins, i, j):
    s_ij =  spin_matrix[i, j]
    s_i1j = spin_matrix[(i+1)%N_spins, j]
    s_ij1 = spin_matrix[i,(j+1)%N_spins]
    s_i0j = spin_matrix[(i-1)%N_spins, j]
    s_ij0 = spin_matrix[i,(j-1)%N_spins]
    dE = 2*s_ij*(s_i1j+s_ij1+s_i0j+s_ij0)
    return dE

def init_config(N_spins):   
    S = 2*np.random.randint(2, size=(N_spins,N_spins))-1
    return S

# compute magnetization
def magnetization(spin_matrix):
    M = np.sum(spin_matrix)
    return M

# compute susceptibility
def susceptibility(m1,m2,N_samples,N_spins,Boltzmann):
    avg_m2 = m2/(N_samples*N_spins**2)
    avg_m1_2 = (m1**2)/((N_samples*N_spins)**2)
    return Boltzmann*(avg_m2 - avg_m1_2)

# approximate critical temperature ased on the observable
def Tc_aprox(X,T):
    ind = X.argsort()[-2:][::-1]
    Tc_sim = np.round(sum(X[ind]*T[ind])/(sum(X[ind])), decimals=2)
    return Tc_sim

def start_simulation(N_spins):
    # get system parameters
    N_spins, N_eq, N_samples, nT, T = parameters(N_spins)
    M = np.zeros(nT)  ; X = np.zeros(nT)

    for n in range(nT):
        print(n)
        Boltzmann=1.0/T[n] ; m1 = 0; m2 = 0
        spin_matrix = init_config(N_spins)

        # Equilibriate
        for i in range(N_eq):         
            spin_matrix = montecarlo(spin_matrix, Boltzmann, N_spins)   

        # MC sample using Metropolis to compute observable
        for i in range(N_samples):
            spin_matrix = montecarlo(spin_matrix, Boltzmann, N_spins)           
            m = magnetization(spin_matrix) 
            m1 = m1 + m
            m2 = m2 + m**2 

        M[n] = m1/(N_samples*N_spins**2)
        X[n] = susceptibility(m1,m2,N_samples,N_spins,Boltzmann)

    return M , X , T

# Simulation parameters
def parameters(N_spins):
    N_spins = N_spins      # Lattice: N_spins x N_spins
    N_eq = 1000            # Steps to equilibriate
    N_samples = 2000       # Steps to sample configurations   
    nT = 60                # Observation points
    T = np.linspace(1.7, 3.0, nT)
    return N_spins,N_eq,N_samples,nT,T

############################# Start simulation #############################
# Run the simulation for 30x30 spins

print('30x30 spins')
M_30 , X_30 , T = start_simulation(N_spins =30)
print('critical temperature (30x30 spins) ~', "{:.2f}".format(Tc_aprox(X_30,T)))

# Run the simulation for 50x50 spins
print('50x50 spins')
M_50 , X_50 , T = start_simulation(N_spins =50)
print('critical temperature (50x50 spins) ~', "{:.2f}".format(Tc_aprox(X_50,T)))

# Save the Magnetization and Susceptibility
np.save('M_30.npy', M_30); np.save('M_50.npy', M_50)
np.save('X_30.npy', X_30); np.save('X_50.npy', X_50)

################################ Plotting ##################################
plt.figure(figsize=(10,6))
plt.scatter(T,M_30, label = '30x30 Spins', linewidth=2)
plt.scatter(T,M_50, label = '50x50 Spins', linewidth=2)
plt.axhline(y=0, color='k',alpha=0.25)
plt.yticks(fontsize=12); plt.xticks(fontsize=12); plt.legend(fontsize=10)
plt.xlabel('$T$',fontsize=14)
plt.ylabel('$<M>$',fontsize=14)
plt.title('Average magnetization vs Temperature',fontsize=14)
plt.grid(alpha=0.5);plt.show()

plt.figure(figsize=(10,6))
plt.scatter(T,abs(M_30), label = '30x30 Spins', linewidth=2)
plt.scatter(T,abs(M_50), label = '50x50 Spins', linewidth=2)
plt.yticks(fontsize=12); plt.xticks(fontsize=12); plt.legend(fontsize=10)
plt.xlabel('$T$',fontsize=14)
plt.ylabel('$|<M>|$',fontsize=14)
plt.title('Absolute average magnetization vs Temperature',fontsize=14)
plt.grid(alpha=0.5);plt.show()

Tc_onsager = 2.27
Tc_sim_30 = Tc_aprox(X_30,T); Tc_sim_50 = Tc_aprox(X_50,T)
plt.figure(figsize=(10,6))
plt.scatter(T,X_30,label = '30x30 Spins', linewidth=2)
plt.scatter(T,X_50,label = '50x50 Spins', linewidth=2)
plt.axvline(x=Tc_onsager, color='g',label = '$T_c \sim$ %.2f (Onsager)' %Tc_onsager, alpha=0.7, linestyle='-.')
plt.axvline(x=Tc_sim_30, label = '$T_c \sim$ %.2f (30x30)' %Tc_sim_30, alpha=0.7, linestyle='-.')
plt.axvline(x=Tc_sim_50, color='orange',label = '$T_c \sim$ %.2f (50x50)' %Tc_sim_50, alpha=0.7, linestyle='-.')
plt.yticks(fontsize=12); plt.xticks(fontsize=12); plt.legend(fontsize=10)
plt.xlabel('$T$',fontsize=14)
plt.ylabel('$\chi_M$',fontsize=14)
plt.title('Susceptibility vs Temperature',fontsize=14)
plt.grid(alpha=0.5);plt.show()
