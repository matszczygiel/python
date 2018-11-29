#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.sparse
import scipy



dt = 0.001    #time step in fs
Nx = 400
Ny = 300
xmin = -20
xmax = 20
ymin = -15
ymax = 15
max_t = 40

Nh = 1

m = 1.
v0 = 0.5
r = 4.
x0  = -10
y0 = 0.
sig = 2.
e = 1.

folder = "1"

kin = 3.81
hbar = 0.6582
mtoh = 0.0864

#potential
def pot(x, y):
#    if x**2+y**2 < r**2:
#        return v0
#    else:
#        return 0.
    return 0.


Vvec = np.vectorize(pot)

E0 = e*v0 
kx0 = np.sqrt(2*m*mtoh*E0 / hbar)
ky0 = 0.
xgrid, dx = np.linspace(xmin, xmax, num=Nx, endpoint=True, retstep=True, dtype=float)
ygrid, dy = np.linspace(ymin, ymax, num=Ny, endpoint=True, retstep=True, dtype=float)
psix = np.exp(- (xgrid - x0)**2 / 2. / sig**2) * np.exp(1j * kx0 * xgrid)
psix *= np.sqrt(1. / sig / np.sqrt(np.pi))
psiy = np.exp(- (ygrid - y0)**2 / 2. / sig**2) * np.exp(1j * ky0 * ygrid)
psiy *= np.sqrt(1. / sig / np.sqrt(np.pi))

psi = np.kron(psix, psiy)

def makeV():
    res = np.array(Nx*Ny)
    for x in range(Nx):
        res[Nx*x:Nx*(x+1)] = Vvec(xgrid[x], ygrid)
    return scipy.sparse.diags(res, [0], shape=(Nx*Ny, Nx*Ny), format='dia')

def makeT():
    res0 = np.full((Nx*Ny), kin *  2. / (m * dx**2))
    res1 = np.full((Nx*Ny-1), - kin / ( m * dx**2))
    Tx = scipy.sparse.diags([res0, res1, res1], [0, 1, -1],
                            shape=(Nx*Ny, Nx*Ny), format='dia' )

    res0 = np.full((Nx*Ny), kin *  2. / (m * dy**2))
    res1 = np.full((Nx*Ny-1), - kin / ( m * dy**2))
    Ty = scipy.sparse.diags([res0, res1, res1], [0, 1, -1], 
                            shape=(Nx*Ny, Nx*Ny), format='dia')

    return scipy.sparse.kron(Tx, Ty, format='dia')

def makeH():
    tx0 = kin *  2. / (m * dx**2)
    tx1 = - kin / ( m * dx**2)
    ty0 = kin *  2. / (m * dy**2)
    ty1 = - kin / ( m * dy**2)
    res0 = np.full((Nx*Ny), tx0 * ty0)
    respot = np.empty(0)
    for x in range(Nx):
        respot = np.append(respot, Vvec(xgrid[x], ygrid))
    res0 += respot
    res1 = np.full((Nx*Ny-1), tx0 * ty1)
    for i in range(1, Nx):
        res1[Ny*i-1] = 0.
    res2 = np.full((Nx*Ny-Ny+1), tx1 * ty1)
    for i in range(0, Nx):
        res2[Ny*i] = 0.
    res3 = np.full((Nx*Ny-Ny), tx1 * ty0)
    res4 = np.full((Nx*Ny-Ny-1), tx1 * ty1)
    for i in range(1, Nx-1):
        res4[Ny*i-1] = 0.

    return scipy.sparse.diags(
                [res0, res1, res1, res2, res2, res3, res3, res4, res4],
                [0, 1, -1, Ny-1, -(Ny-1), Ny, -Ny, Ny+1, -(Ny+1)], 
                shape=(Nx*Ny, Nx*Ny), format='dia')



time = 0.

name = folder + "/e" + str(round(e, 2)) + "psi"
nStr = str(int(time*1000))
nStr=nStr.rjust(5, '0')
#plt.savefig(name+nStr+".png")
psimat = np.mat(np.absolute(psi))
psimat = np.transpose( psimat.reshape((Nx, Ny)) )
psimax = abs(psimat).max()
psimin = 0.
plt.imshow(psimat, interpolation='bilinear',
               origin='lower', extent=[xmin, xmax, ymin, ymax],
               vmax=psimax, vmin=psimin)

plt.savefig(name+nStr+".png")

H = makeH()

def propagate(ps):
    res = ps.copy()
    for n in range(1, Nh+1):
        res = ps.copy() + (-1j * dt ) * H.dot(res)/ (n * hbar)
    return res


while time < max_t:
    time += dt
    time = round(time, 3)
    print("time: " + str(time))    

    psi = propagate(psi)

    plt.clf()
    nStr = str(int(time*1000))
    nStr=nStr.rjust(5, '0')
    psimat = np.mat(np.absolute(psi))
    psimat = np.transpose( psimat.reshape((Nx, Ny)) )
    plt.imshow(psimat, interpolation='bilinear',
                   origin='lower', extent=[xmin, xmax, ymin, ymax],
                   vmax=psimax, vmin=psimin)
    plt.savefig(name+nStr+".png")


