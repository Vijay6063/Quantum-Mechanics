import numpy as np
import matplotlib.pyplot as plt
import math

hbar = 1.0545718e-34
m = 9.109e-31
e = 1.602e-19

a = int(input("Enter well width in nm: ")) * 1e-9
V0 = int(input("Enter potential depth in eV: ")) * e
n = int(input("Enter number of energy states: "))

N = 300
x = np.linspace(-2*a, 2*a, N)
dx = x[1] - x[0]

V = np.zeros(N)

for i in range(N):
    if abs(x[i]) <= a/2:
        V[i] = -V0

H = np.zeros((N-2, N-2))

for i in range(N-2):

    H[i, i] = (hbar**2 / (m*dx**2)) + V[i+1]

    if i > 0:
        H[i, i-1] = -hbar**2 / (2*m*dx**2)

    if i < N-3:
        H[i, i+1] = -hbar**2 / (2*m*dx**2)


energy, wave = np.linalg.eigh(H)

print("\nFinite Potential Well")



plt.figure(figsize=(9, 5))

plt.plot(x*1e9, V/e, 'k', label="Potential")

count = 0

for i in range(len(energy)):

    if energy[i] < 0 and count < n:

        psi = wave[:, i]

     
        psi = psi / math.sqrt(np.sum(psi**2)*dx)

        print("\nEnergy State:", count+1)
        print("Numerical Energy =",
              round(energy[i]/e, 5), "eV")

        plt.plot(x[1:-1]*1e9,
                 psi*0.5e-9 + energy[i]/e,
                 label="State "+str(count+1))

        count += 1

if count == 0:
    print("No bound states found.")

print("\nTotal bound states =", count)

plt.title("Particle in a Finite Potential Well")
plt.xlabel("Position (nm)")
plt.ylabel("Energy (eV)")
plt.legend()
plt.grid()
plt.show()