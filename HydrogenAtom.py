# Computed energy levels and wavefunctions of Hydrogrn Atom using the Finite Difference Method.

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
#  contants
hbar=1
u=1
e=1
k=1

# 1. Computational Domains and Grid, N=Number of grid points,dr=Step size
N = 1000  
a = 0
b = 20 
r = np.linspace( a,b,N)  
dr = r[1] - r[0] 

# 2. Quantum Number(l) and defening the Potential
l=0
def Vpot():
  return (l * (l + 1))*hbar*hbar / (2*u*r**2) - e*e*k /r  # Effective Coulomb potential

# 3. Kinetic Energy Matrix,Potential Matrix and the Hamiltonian Matrix
T=np.zeros((N-2,N-2))
# diagonal elements
for i in range(N-2):
    for j in range(N-2):
        if i==j:
         T[i,j]=hbar**2/(u*dr**2) 
# lower diagonal
    if i>0: 
         j=i-1
         T[i,j]=-hbar**2/(2*u*dr**2 )
# N-3 upper diagonal
    if i<N-3:
         j=i+1
         T[i,j]=-hbar**2/(2*u*dr**2)
print(T)

V=np.zeros((N-2,N-2))
for i in range (N-2):
  for j in range (N-2):
      if i==j:
        V[i,j]=Vpot()[i+1]
print(V)

H= T+ V



# 4. Solving the Eigenvalue Problem
E,Radial_Wavefunction = eigh(H)
print("Lowest three energies:")
print(E[0:3])

# 5. Ground-state radial wavefunction
u0 = Radial_Wavefunction[:, 0]
# Convert discrete normalization to continuum normalization
u0 /= np.sqrt(np.sum(np.abs(u0)**2) * dr)
# Add boundary conditions
u_full = np.zeros(N)
# Fix boundary conditions: u(0) = 0 and u(r_max) = 0 )
u_full[1:-1] = u0

# 6.Plottig of r v/s u(r) with energy having unit Ha
plt.figure(figsize=(8, 5))
plt.plot(r, u_full, label=f"$E_0 = {E[0]:.3f}$ Ha")
plt.title("Hydrogen Atom Ground-State Radial Wavefunction")
plt.xlabel("r(Bohr radius)")
plt.ylabel("u(r)")
plt.legend()
plt.grid(True)
plt.show()