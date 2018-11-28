import numpy as np
import matplotlib.pyplot as plt

dt = 0.1    #time step in fs
N = 7000
xmin = -200
xmax = 200
max_t = 40

m = 0.3
v0 = 0.4
d = 30
x0  = -100
sig = 10
e = 1.5

folder = "p3"

kin = 3.81
hbar = 0.6582
mtoh = 0.0864

#potential
def pot(x):
    if x < 0:
        return 0.
    elif x > d:
        return 0.
    else:
        return v0


V = np.vectorize(pot)

E0 = e*v0 
k0 = np.sqrt(2*m*mtoh*E0 / hbar)
xgrid, dx = np.linspace(xmin, xmax, num=N, endpoint=True, retstep=True, dtype=float)
psi = np.exp(- (xgrid - x0)**2 / 2. / sig**2) * np.exp(1j * k0 * xgrid)
psi *= np.sqrt(1. / sig / np.sqrt(np.pi))

def makeU():
    pottemp = np.exp(- 1j*dt*V(xgrid) / 2. / hbar)
    res = np.empty((N, N), dtype=complex)
    for x in range(N):
        res[x,:] = pottemp[x] * pottemp * np.exp(1j *m*mtoh*(xgrid[x] - xgrid)**2 / 2. / dt)
    res *= np.sqrt(-1j *m*mtoh / (2 * np.pi *dt))*dx
    return res       

U = makeU()

time = 0.

name = folder + "/e" + str(round(e, 2)) + "psi"
nStr = str(int(time*100))
nStr=nStr.rjust(5, '0')
plt.plot(xgrid, np.absolute(psi))
plt.savefig(name+nStr+".png")

bottom, top = plt.ylim()  # return the current ylim

while time < max_t:
    time += dt
    time = round(time, 2)
    print("time: " + str(time))    
    psi = U @ psi

    plt.clf()
    nStr = str(int(time*100))
    nStr=nStr.rjust(5, '0')
    plt.plot(xgrid, np.absolute(psi))
    plt.ylim((bottom, top)) 
    plt.savefig(name+nStr+".png")


