import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

nx = 520
ny = 180
u0 = 0.04
Re = 220

vis = u0 * ny / float(Re * 2)
tau = 3 * vis + 0.5
q = 9

c = np.array([(x, y) for x in [0, -1, 1] for y in [0, -1, 1]])
i1 = np.arange(q)[np.asarray([ci[0]<0  for ci in c])] 
i2 = np.arange(q)[np.asarray([ci[0]==0 for ci in c])] 
i3 = np.arange(q)[np.asarray([ci[0]>0  for ci in c])] 

noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)]

cx = nx / 4.
r = ny / 2.
obstacle = np.fromfunction(lambda x, y: abs(x - cx)+abs(y) < r, (nx, ny))
obstacle[:,0] = np.ones(nx, dtype=bool)
obstacle[:,-1] = np.ones(nx, dtype=bool)


w = [4/9., 1/9., 1/9., 1/9., 1/36., 1/36., 1/9., 1/36., 1/36.]
def equilibrium(rho, mean_v):
    res = np.ones((q, nx, ny))
    cv = np.dot(c, mean_v.transpose(1, 0, 2))
    u2 = mean_v[0]**2 + mean_v[1]**2
    for i in range(q):
        res[i] += 3 * cv[i] + 9 /2. * cv[i]**2 - 1.5 * u2
        res[i] *= w[i] * rho

    return res

vel = np.fromfunction(lambda d,x,y: (1-d)*u0, (2, nx, ny))

feq = equilibrium(1.0, vel)
fin = feq.copy()

max_time = 200001
time = 0

while time < max_time:
    print("time: " + str(time))

    rho = np.sum(fin, axis=0)           
    u = np.dot(c.transpose(), fin.transpose((1,0,2)))/rho
    u[:,0,:] =vel[:,0,:]
    rho[0,:] = ( np.sum(fin[i2,0,:], axis=0) + 2.*np.sum(fin[i1,0,:], axis=0) ) / (1.-u[0,0,:])
    feq = equilibrium(rho,u)
    fin[i3,0,:] = feq[i3,0,:]

    fout = fin -(fin - feq) / tau

    for i in range(q): 
        fout[i, obstacle] = fin[noslip[i], obstacle]
    
    for i in range(q): 
        fin[i,:,:] = np.roll(np.roll(fout[i,:,:], c[i,0], axis=0), c[i,1], axis=1)
    
    if (time%100==0): 
        plt.clf()
        plt.imshow(np.sqrt(u[0]**2 + u[1]**2).transpose(), cmap=cm.Reds)
        plt.savefig("vel."+str(time/100).zfill(4)+".png")

    time += 1