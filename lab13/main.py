#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt


beta2 = 0.5
g = -1
k = 3

def D(x):
    return -beta2 * x**2

def one_step(dz, D_function, g, tvals, dt, psi):
    omegas = np.fft.fftfreq(tvals.shape[-1], dt) * 2 * np.pi
    
    ft = np.fft.fft(psi, norm="ortho")
    ft *= np.exp(-1j * D_function(1j*omegas) * dz / 2.)
    psi_half = np.fft.ifft(ft, norm="ortho")

    ft = np.fft.fft(psi_half * np.exp(-1j * g * np.abs(psi_half)**2 * dz), norm="ortho")
    ft *= np.exp(-1j * D_function(1j*omegas) * dz / 2.)
    return np.fft.ifft(ft, norm="ortho")


def U(tvals, k):
    return np.sqrt(2 * k) / np.cosh(np.sqrt(2 * k) * tvals)



tt, dt = np.linspace(-10, 10, num = 1000, retstep=True)
psi0 = U(tt, k)


dz = 0.01

for i in range(1000):
    print(i)
    plt.clf()
    if i != 0:
        psi0 = one_step(dz, D, g, tt, dt, psi0)

    psiA = np.exp(1j * k * i * dz) * U(tt, k)
    
    plt.ylim((-3, 3))

    plt.plot(tt, psi0.real)
    plt.plot(tt, psi0.imag)

    plt.plot(tt, psiA.real, 'x', markersize=0.6)
    plt.plot(tt, psiA.imag, 'x', markersize=0.6)

    nStr = str(i)
    nStr=nStr.rjust(3, '0')

    plt.savefig("soliton" + nStr + ".png")
