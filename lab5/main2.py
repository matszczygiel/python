import numpy as np
import matplotlib.pyplot as plt

v0 = 0.2
def V(x):
    a = 40
    if abs(x) > a /2:
        return 0.
    else:
        return -v0

Vvec = np.vectorize(V)

plt.clf()

#parameters
m = 0.1
a = -100
b = 100
N = 500

grid, step = np.linspace(a, b, num=N, retstep=True)

kin = 3.81

H = np.zeros((N, N))
for i in range(N):
    H[i, i] = kin *  2. / (m * step**2) + V(grid[i])
    if i != N-1:
        H[i, i+1] = - kin / ( m * step**2)
    if i != 0:
        H[i, i-1] = - kin / ( m * step**2)

lowest = 3
e, v = np.linalg.eigh(H)

plt.plot(grid, Vvec(grid))
plt.hlines(e[0:lowest], a, b)

for i in range(lowest):
    plt.plot(grid, np.ones(N)*e[i] + 0.8 * v0 * v[:,i] / np.amax(abs(v[:,i]))  )

plt.savefig("potential3.png")
