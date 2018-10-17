import numpy

Ngrid = 1000
p = 0.5
M = 4

max_particles = 10000 + 1

x_mid = int(Ngrid / 2)
y_mid = int(Ngrid / 2)

agregate =  numpy.zeros((Ngrid, Ngrid))
agregate[x_mid][y_mid] = -1
agregate[x_mid][y_mid + 1] = 1
agregate[x_mid][y_mid - 1] = 1
agregate[x_mid + 1][y_mid] = 1
agregate[x_mid - 1][y_mid] = 1

def get_radius(agr):
    res = 0
    for x in range(Ngrid):
        for y in range(Ngrid):
            if agr[x][y] == -1:
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

    Ra = get_radius(agregate)
    Rgen = Ra + 5
    Rkill = Rgen + 150
    if  2 * Rkill > Ngrid:
        Rkill = Ngrid / 2 - 1 

    walker = release_walker(Rgen)
    print("released particle: " + str(particle))

    while True:
        w = random_move(walker[0], walker[1])

        if agregate[w[0]][w[1]] > 0:
            radn = numpy.random.random()
            if radn < p:
                if agregate[w[0]][w[1]] == M:
                    agregate[w[0]][w[1]] = -1

                    if agregate[w[0]][w[1] + 1] == 0:
                        agregate[w[0]][w[1] + 1] = 1

                    if agregate[w[0]][w[1] - 1] == 0:
                        agregate[w[0]][w[1] - 1] = 1

                    if agregate[w[0] + 1][w[1]] == 0:
                        agregate[w[0] + 1][w[1]] = 1

                    if agregate[w[0] - 1][w[1]] == 0:
                        agregate[w[0] - 1][w[1]] = 1

                    print("sticked at: " + str(w[0]) + " " + str(w[1]))
                    cir = Circle(w, radius=1,color="black")
                    a.add_patch(cir)

                else: 
                    agregate[w[0]][w[1]] += 1
                    print("killed at sticky point: " + str(w[0]) + " " + str(w[1]))

                break
            else:
                continue


        if ((w[0] - x_mid)**2 + (w[1]  - y_mid)**2) > Rkill**2:
            print("killed at: " + str(w[0]) + " " + str(w[1]))
            break

        walker = w

    if particle % 500 == 0:
        plt.plot()                                         # plot it
        F.set_size_inches((30, 30))            # physical size of the plot
        nStr=str(particle)
        nStr=nStr.rjust(6, '0') 
        plt.savefig('plot' + nStr + '.png')  

