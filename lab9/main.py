#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import *
from sympy.physics.wigner import wigner_3j
import numpy.linalg as npl

De = 0.001174
Re = 10.14
de = 0.24

mvcmtoau = 1 / 5.14220652 * 10**(-3) 

utoau = 1822.88
mLi = 6.015121 * utoau
mRb = 86.909187 * utoau

KtoH = 3.1668114e-6

Rmax = 50
R0 = 9.5

mu = 1/ (1 / mLi + 1 / mRb)
lmaximal = 10

e = 1.0*10**(-6) * KtoH

def pot_int(x):
    return De * ((Re/x)**12 - 2* (Re/x)**6)

def pot_cent(x, ll):
    return ll*(ll+1) / (2* mu * x**2) 

def dipole(x):
    return de * (Re/x)**6

def w_term(lmax, x, F):
    res = np.zeros((len(x), lmax+1, lmax+1), dtype=float)
    for lx in range(lmax+1):
        for ly in range(lmax+1):
            if lx == ly:
                res[:, lx, ly] = pot_cent(x, lx) + pot_int(x)
            else:
                res[:, lx, ly] = - dipole(x) * F * np.sqrt((2*lx+1)*(2*ly+1)) * wigner_3j(lx, ly, 1, 0, 0, 0)**2
    return res

def k_sqrt(x, lmax, e, F):
    return 2 * mu * (np.identity(lmax+1) * e - w_term(lmax, x, F))

def numerov(e, lmax,  F):
    dr = 2 * np.pi / np.sqrt(k_sqrt(np.array([Re]), 0, e, F)[0, 0, 0]) / 50
    r = np.arange(R0, Rmax, dr)
    k2 = k_sqrt(r, lmax, e, F)

    psi = np.zeros((len(r), lmax+1, lmax+1), dtype=float)
    psi[1,:,:] = np.identity(lmax+1) 

    for i in range(1, len(r)-1):
        if i != 1:
            fnm_inv = psi[i-1,:,:] @ npl.inv(psi[i,:,:])
            fn =  (2 * (np.identity(lmax+1)-5*dr**2 * k2[i,:,:]/12) -(np.identity(lmax+1)+dr**2*k2[i-1,:,:]/12) @ fnm_inv )
            fn = npl.inv(np.identity(lmax+1)+ dr**2 * k2[i+1,:,:]/12) @ fn
            psi[i+1,:,:] = fn @ psi[i,:,:] 
        else:
            fn = (2 * (np.identity(lmax+1)-5*dr**2 * k2[i,:,:]/12))
            fn = npl.inv(np.identity(lmax+1)+ dr**2 * k2[i+1,:,:]/12) @ fn
            psi[i+1,:,:] = fn @ psi[i,:,:]

    Flast = psi[-1,:,:] @ npl.inv(psi[-2,:,:])
    kk = np.sqrt(2*e)

    jnm = np.zeros((lmax+1, lmax+1), dtype=float)
    for l in range(lmax+1): jnm[l, l] = kk*r[-2]* spherical_jn(l, kk*r[-2])
    jn = np.zeros((lmax+1, lmax+1), dtype=float)
    for l in range(lmax+1): jn[l, l] = kk*r[-1]* spherical_jn(l, kk*r[-1]) 
    ynm = np.zeros((lmax+1, lmax+1), dtype=float)
    for l in range(lmax+1): ynm[l, l] = kk*r[-2]* spherical_yn(l, kk*r[-2]) 
    yn = np.zeros((lmax+1, lmax+1), dtype=float)
    for l in range(lmax+1): yn[l, l] = kk*r[-1]* spherical_yn(l, kk*r[-1]) 

    K = npl.inv(Flast @ ynm - yn) @ (jn - Flast @ jnm)
    S = npl.inv(np.identity(lmax+1) - 1j*K) @ (np.identity(lmax+1) + 1j*K)
    sigma_el = np.pi / kk**2 * np.abs(1-np.diagonal(S))**2
    sigma_in = np.pi / kk**2 * (1 - np.abs(np.diagonal(S)))**2
    sc_length = 1. / (1j *kk) * (1 - np.diagonal(S)) / (1 + np.diagonal(S)) 
    
    return sigma_el, sigma_in, sc_length
    
Fs = np.linspace(0, 10, num=50, dtype=float) * mvcmtoau
sig_el = np.empty(len(Fs), dtype=object)
sig_in = np.empty(len(Fs), dtype=object)
sc_len = np.empty(len(Fs), dtype=object)
for f in range(len(Fs)):
    print(Fs[f])
    sig_el[f], sig_in[f], sc_len[f] = numerov(e, lmaximal, Fs[f]) 

plt.xscale('log')
plt.plot(Fs, sig_el[:,0], label="elastic")
plt.plot(Fs, sig_in[:,0], label="inelastic")
plt.legend()
plt.savefig("crosssection.png")

plt.xscale('log')
plt.clf()
plt.plot(Fs, sc_len[:,0], label="scatering length")
plt.legend()
plt.savefig("sc_lenght.png")



