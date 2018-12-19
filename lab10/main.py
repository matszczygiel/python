#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npl

gnLi = -0.000447
gnRb = -0.000294
ge = 2.002319

me = 1.0
mp = 1836.15 * me

mhztoau = 1.5198304e-10

ALi = 152.1  * mhztoau
ARb = 1011.9 * mhztoau

ILimult = 3
IRbmult = 6

muB = 1. / 2. / me
muN = 1. / 2. / mp

def get_spin_ops(mult):
    l = (mult - 1) / 2.
    Sz = -np.identity(mult, dtype=float) * l
    for i in range(mult):
        Sz[i, i] += i 

    Sp = np.zeros((mult, mult), dtype=float)
    for i in range(0, mult-1):
        m = -l + i 
        Sp[i+1, i] = np.sqrt(l*(l+1) - m*(m+1))
    
    Sm = np.transpose(Sp)
    return Sz, Sp, Sm


def get_H(A, gn, imult, B):
    sz, sp, sm = get_spin_ops(2)
    iz, ip, im = get_spin_ops(imult)

    Hhfi  = np.kron(iz, sz) + 0.5 * (np.kron(ip, sm) + np.kron(im, sp))
    Hhfi *= A
    H = ge * muB * B * np.kron(np.identity(imult), sz) - gn * muN * B * np.kron(iz, np.identity(2)) + Hhfi
    return H

gausstoau = 1. / 2.35e+9
bfield = np.linspace(0, 2000, num=150, dtype=float)
evalsLi = np.empty((len(bfield), ILimult*2), dtype=float)
evalsRb = np.empty((len(bfield), IRbmult*2), dtype=float)

for i in range(len(bfield)):
    print(i)
    evalsLi[i,:] = npl.eigvalsh(get_H(ALi, gnLi, ILimult, bfield[i] * gausstoau)) / mhztoau
    evalsRb[i,:] = npl.eigvalsh(get_H(ARb, gnRb, IRbmult, bfield[i] * gausstoau)) / mhztoau

plt.clf()
for i in range(ILimult*2):
    plt.plot(bfield, evalsLi[:,i])
plt.savefig("li.png")

plt.clf()
for i in range(IRbmult*2):
    plt.plot(bfield, evalsRb[:,i])
plt.savefig("rb.png")












