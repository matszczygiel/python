import numpy as np
import scipy
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import operator

lx = 300
ly = 200

I = 1
D = 10
E = 0.05
M = 200
R = 100

max_iters = 20001

land = np.zeros((lx, ly))

for x in range(lx):
    for y in range(ly):
        land[x, y] = float(I * y)


def find_counter(land):
    counter = np.zeros((lx, ly))
    for x in range(lx):
        for y in range(ly):
            xi = x
            yi = y
            counter[xi, yi] += 1
            while(True):
                if yi == 0:
                    break

                h_curr = land[xi, yi]
                dx1 = h_curr - land[(xi + 1) % lx, yi]
                dx2 = h_curr - land[(xi - 1) % lx, yi]    
                dy2 = h_curr - land[xi, yi -1]
                dy1 = 0

                if yi == ly - 1:
                    dy1 = -1
                else:
                    dy1 = h_curr - land[xi, yi + 1]

                llist = [dx1, dx2, dy1, dy2]

                index = 0
                val = llist[0]
                equal = [0]
                for i in range(1, 4):
                    if llist[i] > val:
                        index = i
                        val = llist[i]
                        equal = [i]
                    elif llist[i] == val:
                        equal.append(i)

                index = random.choice(equal)

                if index == 0:
                    xi = (xi + 1) % lx
                elif index == 1:
                    xi = (xi - 1) % lx 
                elif index == 2:
                    yi += 1
                else:
                    yi -= 1

                counter[xi, yi] += 1
    
    return counter



def find_rivers(land):
    counter = find_counter(land)
    rivers = []
    for x in range(lx):
        for y in range(ly):
            if counter[x, y] > R:
                rivers.append((x, y))

    return rivers



def plot_rivers(riv, iter):
    plt.clf()
    F = plt.gcf()
    a = plt.gca()
    plt.xlim((0, lx))
    plt.ylim((0, ly))

    for x in riv:
        cir = Circle(x, radius=0.5,color='black')
        a.add_patch(cir)

    plt.plot()                                         
    F.set_size_inches((lx / 50, ly / 50))
    nStr=str(iter)
    nStr=nStr.rjust(6, '0') 
    plt.savefig('rivers' + nStr + '.png')  


def plot_land(lnd, iter):
    plt.clf()
    fig, ax = plt.subplots()
    im = ax.imshow(np.transpose(land), origin='lower')
    plt.colorbar(im)
    nStr=str(iter)
    nStr=nStr.rjust(6, '0') 
    plt.savefig('land' + nStr + '.png')


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

    wet = np.zeros((lx, ly))
    wet[xi, yi] = 1

    while(True):
        if yi == 0:
            break

        p = []
        h_curr = land[xi, yi]
        dx1 = h_curr - land[(xi + 1) % lx, yi]
        dx2 = h_curr - land[(xi - 1) % lx, yi]    
        dy2 = h_curr - land[xi, yi -1]
        dy1 = 0

        if yi == ly - 1:
            dy1 = -1
        else:
            dy1 = h_curr - land[xi, yi + 1]

        if dx1 >= 0:
            p.append(scipy.exp(E * dx1))
        else:
            p.append(0)

        if dx2 >= 0:
            p.append(scipy.exp(E * dx2) + p[0])
        else:
            p.append(p[0])

        if dy1 >= 0:
            p.append(scipy.exp(E * dy1) + p[1])
        else:
            p.append(p[1])

        if dy2 >= 0:
            p.append(scipy.exp(E * dy2) + p[2])
        else:
            p.append(p[2])

        p_tot = p[3]
        for i in range(4):
            p[i] /= p_tot

#        print(dx1, dx2, dy1, dy2)
#        print(p)


        rand = random.random()
        if p[0] > rand:
            xi = (xi + 1) % lx
        elif p[1] > rand:
            xi = (xi - 1) % lx 
        elif p[2] > rand:
            yi += 1
        else:
            yi -= 1

        wet[xi, yi] = 1

    land[wet == 1] -=D

#    print(land)
    land = perform_avalanches(land)       
#    print(land)

    if iter % 1000 == 0:
        print("finding rivers")
        rivers = find_rivers(land)
        print("ploting")
        plot_rivers(rivers, iter)
        plot_land(land, iter)




