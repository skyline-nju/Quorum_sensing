import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    rho0 = 10.
    Lx = Ly = 2.
    N = int(Lx * Ly * rho0)
    x = np.random.rand(N) * Lx
    y = np.random.rand(N) * Ly
    plt.plot(x, y, ".")
    plt.show()
    plt.close()

    count = 0
    count2 = 0.
    for i in range(x.size):
        dx = x[i] - 1
        dy = y[i] - 1
        dd = dx ** 2 + dy ** 2
        if (dd) < 1:
            count += 1
            count2 += 1 - np.sqrt(dd)
    print(count / np.pi)
    print(count2 / np.pi * 3)

