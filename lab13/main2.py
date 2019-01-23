#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz


beta2 = 0.5
g = -1
k = 3

def D(x):
    return -beta2 * x**2

def V(x):
    return 10 * x**2

def one_step(dz, D_function, g, tvals, dt, psi, pot):
    omegas = np.fft.fftfreq(tvals.shape[-1], dt) * 2 * np.pi
    
    ft = np.fft.fft(psi, norm="ortho")
    ft *= np.exp(-1j * D_function(1j*omegas) * dz / 2.)
    psi_half = np.fft.ifft(ft, norm="ortho")

    ft = np.fft.fft(psi_half * np.exp(-1j *( g * np.abs(psi_half)**2)* dz), norm="ortho")
    ft *= np.exp(-1j * D_function(1j*omegas) * dz / 2.)
    return np.fft.ifft(ft, norm="ortho")


def U(tvals, k):
    return np.sqrt(2 * k) / np.cosh(np.sqrt(2 * k) * tvals)

N = 1000
tt, dt = np.linspace(-10, 10, num = N, retstep=True)
#psi0 = np.exp(-4* tt**2 )
psi0 = U(tt, k)

norm = np.sqrt(trapz(psi0**2, tt))

psiAm = U(tt, k)

dz = 0.0001

for i in range(1000):
    print(i)
    plt.clf()
    if i != 0:
        psi0 = one_step(1j*dz, D, g, tt, dt, psi0, V)

    psi0 *= norm / np.sqrt(trapz(np.abs(psi0)**2, tt))
    #psiA = np.exp(1j * k * i * dz) * U(tt, k)
    
    plt.ylim((0, 3))
    #plt.xlim((-5, 5))

    plt.plot(tt, np.abs(psi0))
    #plt.plot(tt, psi0.imag)

    #plt.plot(tt, psiA.real, 'x', markersize=0.6)
    #plt.plot(tt, psiA.imag, 'x', markersize=0.6)
    plt.plot(tt, psiAm, 'x', markersize=0.6)

    #plt.plot(tt, V(tt))

    nStr = str(i)
    nStr=nStr.rjust(3, '0')

    plt.savefig("imt" + nStr + ".png")
