# solved the Schrodinger Equation numerically using the Finite Difference Method.

import numpy as np
import matplotlib.pyplot as plt
# constants
hbar=1
m=1
w=1

#1.defening the Harmonic Potential
def VPot(x):
  return 0.5*m*w**2*x**2

#2. Computationals Domaibs and Grids,N=no of grid points,a=upper limit,b=lower limit
N=100
a=-5
b=5
x=np.linspace(a,b,N)
dx=x[1]-x[0]


#3.  Kinetic energy matrix,Potential energy matrix and Hamiltonian matrix
T=np.zeros((N-2,N-2))
# diagonal elements
for i in range(N-2):
    for j in range(N-2):
        if i==j:
         T[i,j]=hbar**2/(m*dx**2) 
# lower diagonal
    if i>0: 
         T[i,i-1]=-hbar**2/(2*m*dx**2)
# N-3 upper diagonal
    if i<N-3:
         T[i,i+1]=-hbar**2/(2*m*dx**2)
print(T)

V=np.zeros((N-2,N-2))
for i in range (N-2):
  for j in range (N-2):
      if i==j:
       V[i,j]=VPot(x[i+1])
print(V)

H=T+V
print(H)


# 4. Solving the Eigenvalue Problem
E,Wavefunction=np.linalg.eigh(H)
print("Eigen Values:",E[:3])

#5. plotting first three states wave functions of Harmonic oscillator.

plt.figure(figsize=(8, 5))

for n in range(3):
    psi=np.zeros(N)
    psi[1:-1]=Wavefunction[:,n]
    if psi[np.argmax(np.abs(psi))]<0:
        psi=-psi
        # psi=psi/normalisation factor
    psi=psi/np.sqrt(np.sum(np.abs(psi)**2))*dx
    plt.plot(x,psi,label=f"n={n}, E={E[n]:.3f}")

# 6. Plotting x v/s psi
plt.xlabel("x")
plt.ylabel("psi")
plt.title("First Three Harmonic Oscillator Wavefunctions ")
plt.legend()
plt.grid()
plt.show()