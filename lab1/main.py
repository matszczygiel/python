import numpy

Ngrid = 1000
p = 0.8
M = 10

max_particles = 1000

x_mid = int(Ngrid / 2)
y_mid = int(Ngrid / 2)

agregate =  numpy.zeros((Ngrid, Ngrid))
agregate[x_mid][y_mid] = 1

def get_sticky_and_R(agr):
    res = numpy.zeros((Ngrid, Ngrid))
    r = 0
    for x in range(Ngrid):
        for y in range(Ngrid):
            if agr[x][y] == 1:
                if agr[x+1][y] == 0:
                    res[x+1][y] = 1
                if agr[x][y+1] == 0:
                    res[x][y+1] = 1
                if agr[x-1][y] == 0:
                    res[x-1][y] = 1
                if agr[x][y-1] == 0:
                    res[x][y-1] = 1
                rr = (x - x_mid)**2 + (y - y_mid)**2
                if rr > r:
                    r = rr

    return res, numpy.sqrt(r)

def get_radius(agr):
    res = 0
    for x in range(Ngrid):
        for y in range(Ngrid):
            if agr[x][y] == 1:
                r = (x - x_mid)**2 + (y - y_mid)**2
                if r > res:
                    res = r
    return numpy.sqrt(res)
                

import random
def random_move(x, y):
    rand = random.randint(1, 4)
    if rand == 1:
        x += 1
    elif rand == 2:
        y += 1
    elif rand == 3:
        x -= 1
    elif rand == 4:
        y -= 1

    return (x, y)




def release_walker(r):
    a = numpy.random.random()
    w = ((int(r * numpy.cos(a * 2. *numpy.pi) + x_mid), int(r * numpy.sin(a * 2. *numpy.pi) + y_mid)))
    return w


import matplotlib.pyplot as plt
from matplotlib.patches import Circle
plt.clf()
F = plt.gcf()
a = plt.gca()
plt.xlim((0, Ngrid))
plt.ylim((0, Ngrid))

cir = Circle((x_mid, y_mid), radius=1,color='black')
a.add_patch(cir)

for particle in range(max_particles):

    sticky, Ra = get_sticky_and_R(agregate)
    Rgen = Ra + 5
    Rkill = Rgen + 150
    if  2 * Rkill > Ngrid:
        Rkill = Ngrid / 2 - 1 
    print(Ra)
    print(Rgen)
    print(Rkill)

    walker = release_walker(Rgen)
    print("particle released: " + str(particle))

    while True:
        walker_new = random_move(walker[0], walker[1])

        if sticky[walker_new[0]][walker_new[1]] > 0:
            radn = numpy.random.random()
            if radn < p:
                agregate[walker_new[0]][walker_new[1]] = 1
                print("sticked at: " + str(walker_new[0]) + " " + str(walker_new[1]))
                cir = Circle(walker_new, radius=1,color="black")
                a.add_patch(cir)
                break
            else:
                continue


        if ((walker_new[0] - x_mid)**2 + (walker_new[1]  - y_mid)**2) > Rkill**2:
            print("killed at: " + str(walker_new[0]) + " " + str(walker_new[1]))
            break

        walker = walker_new

    if particle % 50 == 0:
        plt.plot()                                         # plot it
        F.set_size_inches((30, 30))            # physical size of the plot
        nStr=str(particle)
        nStr=nStr.rjust(6, '0') 
        plt.savefig('plo' + nStr + '.png')  

