#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

sig = 1.

N = int(2**10)

def gaussian(x, sigma):
    res = np.exp(-0.5 * (x / sigma)**2 )
    res /= np.sqrt(sigma) * np.pi**0.25
    return res

xv, dx = np.linspace(-40, 40, num=N, dtype=float, retstep=True)
gv = gaussian(xv, sig)

gvf = np.fft.fft(gv, norm="ortho")
xvf = np.fft.fftfreq(xv.shape[-1], dx) * 2 * np.pi

phasef = np.exp(-1j * xv[0] * xvf)

gvf *= phasef

plt.clf()
plt.plot(xv, gv)
plt.savefig("gv.png")

plt.clf()
plt.plot(xvf, gvf.real, 'o')
plt.plot(xvf, gaussian(xvf, 1./sig), 'x')
plt.savefig("gvf.png")



