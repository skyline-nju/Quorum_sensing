import numpy as np
import matplotlib.pyplot as plt


def get_slope(p1, p2):
    return (p1[1] - p2[1]) / (p1[0] - p2[0])


if __name__ == "__main__":
    # two-phase coexisting
    rc = [[7.5, 10], [6.5, 14.05]]

    p_g = [[3.8842, 1.1965], [3.4236, 1.3114]]
    p_l = [[10.5404, 21.6501], [7.1511, 21.9117]]

    # k0 = (r1[0][1] - r2[0][1])/(r1[0][0] - r2[0][0])
    k0 = get_slope(p_g[0], p_l[0])
    k1 = -1 / k0


    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    for i in range(len(rc)):
        line, = ax.plot(rc[i][0], rc[i][1], "s")
        ax.plot(p_g[i][0], p_g[i][1], "v", fillstyle="none", c=line.get_c())
        ax.plot(p_l[i][0], p_l[i][1], "^", fillstyle="none", c=line.get_c())
        ax.plot([p_g[i][0], p_l[i][0]], [p_g[i][1], p_l[i][1]], linestyle="dotted", c=line.get_c())


    # three-phase coexisting
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
    x2 = 3.5
    y2 = (x2-x1) * k1 + y1
    
    x3 = 10.5
    y3 = (x3-x1) * k1 + y1

    ax.plot(x1, y1, ">")
    ax.plot(x2, y2, ">")
    ax.plot(x3, y3, ">")

    ax.plot(4.5, 6, "p")
    ax.plot(16.6, 13.6, "p")


    x = 4.9
    k = get_slope(p_g[1], p_l[1])
    y = 14.05 + k * (x-6.5)
    ax.plot(x, y, "^")
    print(x, y, k)

    ax.axis("scaled")
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 25)
    print(x2, y2)
    print(x3, y3)
    plt.show()
    plt.close()
        
