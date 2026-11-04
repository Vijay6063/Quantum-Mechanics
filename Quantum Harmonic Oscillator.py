import numpy as np
import matplotlib.pyplot as plt

l = float(input("Enter length limit (l): "))
n = int(input("Enter number of grid points (n): "))
m = float(input("Enter mass of particle (m): "))
hbar = float(input("Enter reduced Planck's constant (hbar): "))
k = float(input("Enter force constant (k): "))

x = np.linspace(-l, l, n)
dx = x[1] - x[0]
V = 0.5 * k * x**2

main_diag = -2 * np.ones(n - 2)
off_diag = 1 * np.ones(n - 3)

K_E = np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)

H = -((hbar**2) / (2 * m * dx**2)) * K_E + np.diag(V[1:-1])

energies, wavefunctions = np.linalg.eigh(H)

for i in range(3):
    psi = wavefunctions[:, i]
    psi = psi / np.sqrt(np.sum(psi**2) * dx)

    psi_x = np.zeros(n)
    psi_x[1:-1] = psi

    psi_y = (psi_x * 4.0) + energies[i]

    plt.plot(x, psi_y, label=f'State n={i} (E={energies[i]:.2f})')

plt.title('Quantum Harmonic Oscillator & Energy Levels')
plt.xlabel('Position (x)')
plt.ylabel('Energy / Wavefunction Amplitude')
plt.plot(x, V, color='black', label='V(x) = 0.5 * k * x^2')
plt.grid(True)
plt.legend()
plt.show()