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
stotvals = np.arange(-4.5, 4.51, step=0.5)
scount = np.zeros(len(stotvals), dtype=int)
evals = np.empty((len(stotvals)), dtype=object)

SzLi = np.kron(get_spin_ops(ILimult)[0], np.identity(2)) + np.kron(np.identity(ILimult), get_spin_ops(2)[0])
SzRb = np.kron(get_spin_ops(IRbmult)[0], np.identity(2)) + np.kron(np.identity(IRbmult), get_spin_ops(2)[0])
Sztot = np.kron(SzLi, np.identity(IRbmult*2)) + np.kron(np.identity(ILimult*2), SzRb)
Sdiag = np.diagonal(Sztot)

for i in range(1, len(bfield)):
    print(i)
    H = np.kron(get_H(ALi, gnLi, ILimult, bfield[i] * gausstoau), np.identity(IRbmult*2))
    H += np.kron(np.identity(ILimult*2), get_H(ARb, gnRb, IRbmult, bfield[i] * gausstoau))
    for s in range(len(stotvals)):
        egvs, vecs = npl.eigvalsh(hblock) 
        egvs /= mhztoau
        spins = np.diagonal(vecs.transpose() @ Sztot @ vecs) 
            if i == 0:
                evals[s] = np.zeros((len(bfield), len(egvs)), dtype=float)
                scount[s] = len(egvs)
            evals[s][i, :] = egvs
        else:
            if i == 0:
                evals[s] = np.zeros((len(bfield), 0), dtype=float)
                scount[s] = 0

print(scount)
for s in range(len(stotvals)):
    plt.clf()

    for e in range(scount[s]):
        plt.plot(bfield, evals[s][:, e])
        print(s, e, evals[s][:, e])
    plt.savefig("plot" + str(stotvals[s]) + ".png")
#
#ev, vecs = npl.eigh(H)
#print(sv)
#print(ev / mhztoau)
#
#for s in np.arange(-4.5, 4.51, step=0.5):
#    plt.clf()
#    evs = evals[:,sv == s]
#    for e in range(len(evs[0,:])):
#        plt.plot(bfield, evs[:,e])
#        plt.savefig("plot" + str(s) + ".png")
#


#plt.clf()
#for i in range(ILimult*IRbmult*4):
#    plt.plot(bfield, evals[:,i])












