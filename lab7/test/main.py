import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.linalg

dt = 0.01    #time step in fs
N = 400
xmin = -20
xmax = 20
max_t = 40

Nh = 70

m = 1
x0  = -10.
sig = 2. 
E0 = 0.5

kin = 3.81
hbar = 0.6582
mtoh = 0.0864

k0 = np.sqrt(2*m*mtoh*E0 / hbar)
xgrid, dx = np.linspace(xmin, xmax, num=N, endpoint=True, retstep=True, dtype=float)
psi = np.exp(- (xgrid - x0)**2 / 2. / sig**2) * np.exp(1j * k0 * xgrid)
psi *= np.sqrt(1. / sig / np.sqrt(np.pi))


H = np.zeros((N, N))
for i in range(N):
    H[i, i] = kin *  2. / (m * dx**2)
    if i != N-1:
        H[i, i+1] = - kin / ( m * dx**2)
    if i != 0:
        H[i, i-1] = - kin / ( m * dx**2)



def propagate(ps):
    tmp = ps.copy()
    scal = -1.j * dt / hbar
    mult = 1. + 0.j
    res = tmp.copy()
    for n in range(1, Nh+1):
        mult *= scal / n
        tmp = H @ tmp
        res += tmp * mult
    return res


time = 0.

name = "e" + str(round(E0, 2)) + "psi"
nStr = str(int(time*100))
nStr=nStr.rjust(5, '0')
plt.plot(xgrid, np.absolute(psi))
plt.savefig(name+nStr+".png")

bottom, top = plt.ylim()  # return the current ylim


while time < max_t:
    time += dt
    time = round(time, 2)
    print("time: " + str(time))    
    
    psi = propagate(psi)


    plt.clf()
    nStr = str(int(time*100))
    nStr=nStr.rjust(5, '0')
    plt.plot(xgrid, np.absolute(psi))
    plt.ylim((bottom, top)) 
    plt.savefig(name+nStr+".png")


