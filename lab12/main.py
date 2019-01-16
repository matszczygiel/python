#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt


beta1 = 1
beta2 = 1

def D(x):
    return -beta1 * x - 1j* beta2 / 2. * x**2

def d_step(dz, D_function, tvals, dt, psi):
    ft = np.fft.fft(psi, norm="ortho")
    omegas = np.fft.fftfreq(tvals.shape[-1], dt) * 2 * np.pi
    ft *= np.exp(D_function(1j*omegas) * dz)
    return np.fft.ifft(ft, norm="ortho")

t0 = 1
A0 = 1

tt, dt = np.linspace(-10, 10, num = 1000, retstep=True)
A0tab = A0 * np.exp(- 0.5 * (tt / t0)**2)

plt.clf()
plt.plot(tt, A0tab)
dz = 1
Aprop = d_step(dz, D, tt, dt, A0tab)
plt.plot(tt, Aprop.real)
plt.plot(tt, Aprop.imag)

def A_analit(tvals, z):
    return A0 * t0 / np.sqrt(t0**2 - 1j * beta2 * z) * np.exp(-0.5 * (tt - beta1 * z)**2 / (t0**2 - 1j* beta2 * z))

AAa = A_analit(tt, dz)
plt.plot(tt, AAa.real, 'x', markersize=0.6)
plt.plot(tt, AAa.imag, 'x', markersize=0.6)


plt.savefig("gauss.png")
