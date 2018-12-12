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

mu = 1/ (1 / m40 + 1 / m39)
lmax = 15

ens = np.linspace(-6, 0, num=300, dtype=float)
Ens = 10**ens
sigtot = np.zeros(len(ens))
sig = np.zeros((len(ens), lmax+1))

for i in range(len(Ens)):
    E = Ens[i]
    e = E * KtoH
    kk = np.sqrt(2*mu*e)
    
    def pot_int(x):
        return De * ((Re/x)**12 - 2* (Re/x)**6)
    
    def pot(x, ll):
        return ll*(ll+1) / (2* mu * x**2) + pot_int(x) 
    
    def k_sqrt(x, ll):
        return 2 * mu * (e - pot(x, ll))
   
    def numerov(l):
        dr = 2 * np.pi / np.sqrt(k_sqrt(Re, l))/ 50
        r = np.arange(R0, Rmax, dr)
        k2 = k_sqrt(r, l)

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
        jnm = kk*r[-2] * spherical_jn(l, kk*r[-2]) 
        jn  = kk*r[-1] * spherical_jn(l, kk*r[-1]) 
        ynm = kk*r[-2] * spherical_yn(l, kk*r[-2]) 
        yn  = kk*r[-1] * spherical_yn(l, kk*r[-1])
        K = (F*jnm  - jn) / (yn - F*ynm)
        s = (1+1j*K) / (1-1j*K)
        sigma = 4*np.pi / kk**2 * (2*l+1) * np.abs(1-s)**2
        return r, psi, sigma
    
    
    
    r = np.empty(lmax+1, dtype=object)
    psi = np.empty(lmax+1, dtype=object)
    
    for l in range(lmax+1):
        r[l], psi[l], sig[i,l] = numerov(l)
    
    sigtot[i] = np.sum(sig[i])

plt.xscale('log')
plt.yscale('log')
plt.ylim(10**(-2), 10**7)
plt.plot(Ens, sigtot, label="total")
for l in range(0,lmax+1):
    ls = str(l)
    plt.plot(Ens, sig[:,l], label="l = "+ls)
plt.legend()
plt.savefig("crosssection.png")


