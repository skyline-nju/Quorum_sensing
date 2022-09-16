import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    rc = [[7.5, 10], [6, 10]]

    r1 = [[3.8842, 1.1965], [3.51, 1.226]]
    r2 = [[10.5404, 21.6501], [7.5327, 22.0631]]

    k0 = (r1[0][1] - r2[0][1])/(r1[0][0] - r2[0][0])
    k1 = -1 / k0


    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    for i in range(len(rc)):
        line, = ax.plot(rc[i][0], rc[i][1], "s")
        ax.plot(r1[i][0], r1[i][1], "o", fillstyle="none", c=line.get_c())
        ax.plot(r2[i][0], r2[i][1], "o", fillstyle="none", c=line.get_c())
        ax.plot([r1[i][0], r2[i][0]], [r1[i][1], r2[i][1]], linestyle="dotted", c=line.get_c())


    rc = [[12.5, 10]]
    r1 = [[4.454, 1.197]]
    r2 = [[13.245, 21.375]]
    r3 = [[21.4, 2.88]]
    for i in range(len(rc)):
        line, = ax.plot(rc[i][0], rc[i][1], "s")
        ax.plot(r1[i][0], r1[i][1], "o", fillstyle="none", c= line.get_c())
        ax.plot(r2[i][0], r2[i][1], "o", fillstyle="none", c= line.get_c())
        ax.plot(r3[i][0], r3[i][1], "o", fillstyle="none", c= line.get_c())


    x1, y1 = 8.625, 13.36
    x2 = 6.5
    y2 = (x2-x1) * k1 + y1
    
    x3 = 10.5
    y3 = (x3-x1) * k1 + y1

    ax.plot(x1, y1, ">")
    ax.plot(x2, y2, ">")
    ax.plot(x3, y3, ">")

    ax.plot(4.5, 6, "p")
    ax.axis("scaled")
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 25)
    print(x2, y2)
    print(x3, y3)

    plt.show()
    plt.close()
        
