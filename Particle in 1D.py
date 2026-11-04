
import numpy as np
import matplotlib.pyplot as plt
import math


hbar = 1.0545718e-34
m = 9.109e-31
e = 1.602e-19


L = int(input("Enter box length in nm: ")) * 1e-9
n = int(input("Enter number of energy states: "))


N = 200
x = np.linspace(0, L, N)
dx = x[1] - x[0]

H = np.zeros((N-2, N-2))

for i in range(N-2):
    H[i, i] = 2

    if i > 0:
        H[i, i-1] = -1

    if i < N-3:
        H[i, i+1] = -1

H = H * (hbar**2 / (2*m*dx**2))


energy, wave = np.linalg.eigh(H)

print("\nEnergy Eigenvalues")

for i in range(n):

    E_num = energy[i] / e

    E_exact = ((i+1)**2 * math.pi**2 * hbar**2) / (2*m*L**2*e)

    print("\nEnergy State:", i+1)
    print("Numerical Energy =", round(E_num, 5), "eV")
    print("Analytical Energy =", round(E_exact, 5), "eV")


plt.figure(figsize=(8, 5))

for i in range(n):

    psi = wave[:, i]

    psi = psi / math.sqrt(np.sum(psi**2)*dx)

    plt.plot(x[1:-1]*1e9, psi, label="State "+str(i+1))

plt.title("Particle in a Box")
plt.xlabel("Position (nm)")
plt.ylabel("Wave Function")
plt.legend()
plt.grid()
plt.show()

