import numpy as np
import matplotlib.pyplot as plt

v0 = 0.2
def V(x):
    a = 40. 
    res = (4. * v0 / a**2) *  x**2 - v0
    if res > 0.:
        return 0.
    else:
        return res

Vvec = np.vectorize(V)

plt.clf()

#parameters
m = 0.1
a = -100
b = 100
N = 501 #must be odd

grid, step = np.linspace(a, b, num=N, retstep=True)

kin = 3.81

kr = np.arange((N + 1)/2)
T = 4 * kin / m * (np.pi / step / N * kr )**2

H = np.zeros((N, N))
for x in range(N):
    for y in range(N):
        H[x, y] = 2. / N * np.sum(T * np.cos(2 * np.pi * (x-y) * kr / N))
        if x == y:
            H[x, x] += V(grid[x])

lowest = 4
e, v = np.linalg.eigh(H)

plt.plot(grid, Vvec(grid))
plt.hlines(e[0:lowest], a, b)

for i in range(lowest):
    plt.plot(grid, np.ones(N)*e[i] + 0.8 * v0 * v[:,i] / np.amax(abs(v[:,i]))  )

plt.savefig("potential4.png")
