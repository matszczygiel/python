#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import *

De = 0.0011141
Re = 10.98

utoau = 1822.88
m39 = 38.963707 * utoau
m40 = 39.693999 * utoau
m41 = 40.961825 * utoau

KtoH = 3.1668114e-6

Rmax = 50
R0 = 9.5

mu = 1/ (1 / m39 + 1 / m39)
lmax = 6

E = 1.e-2
e = E * KtoH
kk = np.sqrt(2*e)

lamb = np.linspace(0.95, 1.05, num=1000)
sigtot = np.zeros(len(lamb))

for i in range(len(lamb)):
    def pot_int(x):
        return lamb[i]*De * ((Re/x)**12 - 2* (Re/x)**6)

    def pot(x, ll):
        return ll*(ll+1) / (2* mu * x**2) + pot_int(x) 

    def k(x, ll):
        return np.sqrt(2 * mu * (e - pot(x, ll)) + 0j)


    def numerov(l):
        dr = 2 * np.pi / k(Re, l).real / 50
        r = np.arange(R0, Rmax, dr)
        k2 = (k(r, l)**2).real

        psi = np.zeros(len(r), dtype=float)
        psi[1] = 1.

        for i in range(1, len(r)-1):
            if i != 1:
                fnm = psi[i] / psi[i-1]
                fn = (2 * (1-5*dr**2 * k2[i]/12) -(1+dr**2*k2[i-1]/12)/fnm )/ (1+ dr**2 * k2[i+1]/12)
                psi[i+1] = psi[i] * fn
            else:
                fn = 2 * (1-5*dr**2 * k2[i]/12)/ (1+ dr**2 * k2[i+1]/12)
                psi[i+1] = psi[i] * fn

        F = psi[-1] / psi[-2]
        jnm = spherical_jn(l, kk*r[-2]) 
        jn = spherical_jn(l, kk*r[-1]) 
        ynm = spherical_yn(l, kk*r[-2]) 
        yn = spherical_yn(l, kk*r[-1])
        K = (F*jnm  - jn) / (yn - F*yn)
        delta = np.arctan(-K)
        sigma = np.pi / (2*e) * (2*l+1) * np.sin(delta)**2
        return r, psi, sigma



    r = np.empty(lmax+1, dtype=object)
    psi = np.empty(lmax+1, dtype=object)
    sig = np.empty(lmax+1)
    for l in range(lmax+1):
        r[l], psi[l], sig[l] = numerov(l)

    sigtot[i] = np.sum(sig)

plt.plot(lamb, sigtot)
#for l in range(1):
#    plt.plot(ens, sig[:,l])
plt.savefig("cross10mk.png")

