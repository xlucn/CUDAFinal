import matplotlib.pyplot as plt
import numpy as np

def ploteach(datafile):
    with open(datafile) as f:
        u = [line.split() for line in f]

    U = np.array(u)
    ly, lx = U.shape

    x = np.linspace(0, 100, lx)
    y = np.linspace(0, 100, ly)
    X, Y = np.meshgrid(x, y)

    plt.contour(X, Y, U)
    plt.title(datafile)
    plt.colorbar()

plt.figure(figsize=(12, 8))
plt.subplot(221)
ploteach("u0.txt")
plt.subplot(222)
ploteach("exact.txt")
plt.subplot(223)
ploteach("diff.txt")
plt.savefig("plot.png")
plt.close()
