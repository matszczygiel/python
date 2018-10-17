import numpy as np 
import random

lx = 10
ly = 10

I = 1
D = 10
E = 0.05
M = 5
R = 20

max_iters = 5

land = np.zeros((lx, ly))
counter = np.zeros((lx, ly))

for x in range(lx):
    for y in range(ly):
        land[x, y] = float(I * y)

def perform_avalanches(l):
    avalanches = False
    while(not avalanches):
        land_tmp = l.copy()

        difference = l - np.roll(l, -1, axis=0)
        land_tmp[difference > M] = land_tmp[difference > M] - difference[difference > M] /4.

        difference = l - np.roll(l, 1, axis=0)
        land_tmp[difference > M] = land_tmp[difference > M] - difference[difference > M] /4.

        difference = l - np.roll(l, -1, axis=1)
        difference[:, -1] = np.zeros(lx)
        land_tmp[difference > M] = land_tmp[difference > M] - difference[difference > M] /4.

        difference = l - np.roll(l, 1, axis=1)
        land_tmp[difference > M] = land_tmp[difference > M] - difference[difference > M] /4.

        avalanches = np.array_equal(land_tmp, l)
        l = land_tmp

    return l
        
for iter in range(max_iters):

    print("Iteration: " + str(iter))

    xi = random.randint(0, lx-1)
    yi = random.randint(0, ly-1)

    print("drop at: " + str(xi) + " " + str(yi))

    wet = [(xi, yi)]

    while(True):
        counter[xi][yi] += 1
        p = []
        h_curr = land[xi, yi]
        dx1 = h_curr - land[(xi + 1) % lx, yi]
        dx2 = h_curr - land[(xi - 1) % lx, yi]    
        dy2 = h_curr - land[xi][yi -1]
        dy1 = 0

        if yi == 0:
            break
        elif yi == ly - 1:
            dy1 = 0
        else:
            dy1 = h_curr - land[xi, yi + 1]
        if dx1 > 0:
            p.append(E * dx1)
        else:
            p.append(0)

        if dx2 > 0:
            p.append(E * dx2)
        else:
            p.append(0)
        if dy1 > 0:
            p.append(E * dy1)
        else:
            p.append(0)
        if dy2 > 0:
            p.append(E * dy2)
        else:
            p.append(0)

        p_tot = np.sum(p)
        for i in range(4):
            p[i] /= p_tot

        rand = random.random()
        if p[0] > rand:
            xi = (xi + 1) % lx
        elif np.sum(p[0:2]) > rand:
            xi = (xi - 1) % lx 
        elif np.sum(p[0:3]) > rand:
            yi += 1
        else:
            yi -= 1

        wet.append((xi, yi))

    for w in wet:
        land[w] -= D

    print(land)
    land = perform_avalanches(land)       
    print(land)
   

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
plt.clf()
F = plt.gcf()
a = plt.gca()
plt.xlim((0, lx))
plt.ylim((0, ly))

for x in range(lx):
    for y in range(ly):
        if counter[x][y] > R:
            cir = Circle((x, y), radius=0.5,color='black')
            a.add_patch(cir)

plt.plot()                                         
F.set_size_inches((10, 10))
plt.savefig('rivers.png')  

plt.clf()
fig, ax = plt.subplots()
im = ax.imshow(np.transpose(land), origin='lower')
plt.savefig('land.png')
