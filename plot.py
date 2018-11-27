import matplotlib.pyplot as plt
import numpy as np

def ploteach(ax, datafile):
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
ax1 = plt.subplot(221)
ploteach(ax1, "u0.txt")
ax2 = plt.subplot(222)
ploteach(ax2, "exact.txt")
ax3 = plt.subplot(223)
ploteach(ax3, "diff.txt")
plt.savefig("plot.png")
plt.close()
