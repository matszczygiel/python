#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

sig = 1.
gamma = 1.
omega0 = 10.
delt = 2.

N = int(2**10)



def compute_fft(xvals, dx, fvals):
    ff = np.fft.fft(fvals, norm="ortho")
    freqf = np.fft.fftfreq(xvals.shape[-1], dx) * 2 * np.pi
    ff *=  np.exp(-1j * xvals[0] * freqf) 
    return freqf, ff 


def gaussian(x, sigma):
    res = np.exp(-0.5 * (x / sigma)**2 )
    res /= np.sqrt(sigma) * np.pi**0.25
    return res

def exponent(x ,gamma):
    return np.exp(-gamma * np.abs(x))

def lorentzian(x, gamma):
    return 2 * gamma / (gamma**2 + x**2)/ np.sqrt(2 * np.pi)

def turn_cos(x, delt, omega0):
    res = np.zeros(x.shape, dtype=float)
    for i in range(len(x)):
        if x[i] > -0.5*delt and x[i] < 0.5*delt:
            res[i] = np.cos(omega0 * x[i])
    return res

def turn_cos_FT(x, delt, omega0):
    res = np.sin(omega0 - x) / (omega0 - x) + np.sin(omega0 + x) / (omega0 + x)
    res *= delt / 2. / np.sqrt(2 * np.pi)
    return res


xv, dx = np.linspace(-40, 40, num=N, dtype=float, retstep=True)
gv = gaussian(xv, sig)
ev = exponent(xv, gamma)
sv = turn_cos(xv, delt, omega0)

xgvf, gvf = compute_fft(xv, dx, gv) 
xevf, evf = compute_fft(xv, dx, ev) 
xsvf, svf = compute_fft(xv, dx, sv) 

mksize = 3

plt.clf()
plt.plot(xv, gv)
plt.savefig("gv.png")

plt.clf()
plt.plot(xgvf, gvf.real, 'o', markersize=mksize)
plt.plot(xgvf, gaussian(xgvf, 1./sig), 'x',markersize=mksize*1.5)
plt.savefig("gvf.png")

plt.clf()
plt.plot(xv, ev)
plt.savefig("ev.png")

plt.clf()
plt.plot(xevf, evf.real, 'o', markersize=mksize)
plt.plot(xevf, lorentzian(xevf, gamma), 'x',markersize=mksize*1.5)
plt.savefig("evf.png")

plt.clf()
plt.plot(xv, sv)
plt.savefig("sv.png")

plt.clf()
plt.plot(xsvf, svf.real, 'o', markersize=mksize)
plt.plot(xsvf, turn_cos_FT(xsvf, delt, omega0), 'x',markersize=mksize*1.5)
plt.savefig("svf.png")


