#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.sparse
import scipy
import scipy.sparse.linalg
import scipy.linalg

dt = 0.1    #time step in fs
Nx = 160
Ny = 160
xmin = -40
xmax = 40
ymin = -40
ymax = 40
max_t = 40

Nh = 100

time_mult = 10

m = 1.
v0 = 0.5
r = 4.
x0  = -10.
y0 = 0.
sig = 4.
e = 1.

folder = "3"

kin = 3.81
hbar = 0.6582
mtoh = 0.0864

#potential
def pot(x, y):
    if x**2+y**2 < r**2:
        return v0
    else:
        return 0.
#    return 0.


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

def makeH():
    t0 = kin *  (2. / m) * (1. / dx**2 + 1. / dy**2)
    tNx = - kin / ( m * dx**2)
    t1 = - kin / ( m * dy**2)
    res0 = np.full((Nx*Ny), t0)
    respot = np.empty(0)
    for x in range(Nx):
        respot = np.append(respot, Vvec(xgrid[x], ygrid))
    res0 += respot
    res1 = np.full((Nx*Ny-1), t1)
    resNx = np.full((Nx*Ny - Nx), tNx)
    return scipy.sparse.diags(
                [res0, res1, res1, resNx, resNx],
                [0, 1, -1, Nx, -Nx], 
                shape=(Nx*Ny, Nx*Ny), format='dia')

    
name = folder + "/e" + str(round(e, 2)) + "psi"
psimin = 0.
psimax = np.absolute(psi).max()
def picture(ps, suffix):
    plt.clf()
    psimat = np.mat(np.absolute(ps))
    psimat = np.transpose( psimat.reshape((Nx, Ny)) )
    plt.imshow(psimat, interpolation='bilinear',
               origin='lower', extent=[xmin, xmax, ymin, ymax],
               vmax=psimax, vmin=psimin)

    plt.savefig(name+suffix+".png")


time = 0.

nStr = str(int(time*time_mult))
nStr=nStr.rjust(5, '0')
picture(psi, nStr)

H = makeH()

def propagate(ps):
#    tmp = ps.copy()
    scal = -1.j * dt / hbar
#    mult = 1. + 0.j
#    res = tmp.copy()
#    for n in range(1, Nh+1):
#        mult *= scal / n
#        tmp = H.dot(tmp * mult)
#        res += tmp        
#    return res
    return scipy.sparse.linalg.expm_multiply(scal * H, ps)

while time < max_t:
    time += dt
    time = round(time, 4)
    print("time: " + str(time))    

    psi = propagate(psi)

    nStr = str(int(time*time_mult))
    nStr=nStr.rjust(5, '0')
    picture(psi, nStr)

