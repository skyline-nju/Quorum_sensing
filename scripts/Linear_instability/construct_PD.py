import numpy as np
import matplotlib.pyplot as plt

# from profile import plot_time_ave_density_profile_SP


def get_slope(p1, p2):
    return (p1[1] - p2[1]) / (p1[0] - p2[0])

def plot_tri(ax, edge, center, c, linestyle="dashed"):
    for i in range(len(edge)):
        ax.plot([edge[i-1][0], edge[i][0]], [edge[i-1][1], edge[i][1]], c=c, linestyle=linestyle)
        ax.plot([edge[i][0], center[0]], [edge[i][1], center[1]], c=c, linestyle="dotted")

def plot_bi(ax, edge, c, linestyle="dotted"):
    ax.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], c=c, linestyle=linestyle)


def plot_tie_lines(ax, centers, edges, marker="o", c=None):
    for i, center in enumerate(centers):
        x = [j[0] for j in edges[i]]
        y = [j[1] for j in edges[i]]
        line, = ax.plot(center[0], center[1], marker, c=c)
        ax.plot(x, y, marker, fillstyle="none", c=line.get_c())
        plot_bi(ax, edges[i], line.get_c())


def fill_color(edge, ax, c="tab:blue"):
    n = len(edge)
    x, y = [], []
    for p in edge:
        x.append(p[0][0])
        y.append(p[0][1])
    for p in edge[::-1]:
        x.append(p[1][0])
        y.append(p[1][1])
    ax.fill(x, y, c=c, alpha=0.5)


def binodal(ax):
    # LB + DC + G
    center_3p = [15/10, 10/10]
    edge_3p = [[13.215/10, 21.3376/10], [4.4581/10, 1.1997/10], [20.6741/10, 2.8125/10]]
    x = [i[0] for i in edge_3p]
    y = [i[1] for i in edge_3p]
    line, = ax.plot(center_3p[0], center_3p[1], ">")
    ax.plot(x, y, ">", fillstyle="none", c=line.get_c())
    plot_tri(ax, edge_3p, center_3p, line.get_c())

    x = [i[0] for i in edge_3p]
    y = [i[1] for i in edge_3p]
    ax.fill(x, y, c="tab:pink", alpha=0.3)

    # DC + G
    center_G_DC = [[9.5/10, 13.075/10], [8.625/10, 13.36/10], [7.5/10, 13.725/10], [6.5/10, 14.05/10]]
    edge_G_DC = [[[4.1251/10, 1.2010/10], [11.8641/10, 21.5710/10]],
                 [[3.9359/10, 1.1988/10], [10.4761/10, 21.9334/10]],
                 [[3.6668/10, 1.2039/10], [8.644/10, 22.1592/10]],
                 [[3.404/10, 1.3104/10], [7.1505/10, 21.9526/10]]
                ]
    plot_tie_lines(ax, center_G_DC, edge_G_DC, "s")
    fill_color([[edge_3p[1], edge_3p[0]]] + edge_G_DC, ax, c="tab:orange")

    # G + LA
    center_G_LA = [
        [12.5/10, 0],
        [12.5/10, 0.5/10],
        [12.5/10, 1/10],
        [12.5/10, 2/10],
        # [12.5, 3]
    ]
    edge_G_LA = [
        [[4.1229/10, 0], [22.1689/10, 0]],
        [[4.1984/10, 0.2967/10], [21.8290/10, 0.7313/10]],
        [[4.2891/10, 0.5983/10], [21.5045/10, 1.4470/10]],
        [[4.4534/10, 1.2062/10], [20.6741/10, 2.8125/10]],
        # [[4.6308, 1.7986], [19.6164, 4.1027]]
    ]
    plot_tie_lines(ax, center_G_LA, edge_G_LA, c="tab:grey")
    fill_color(edge_G_LA, ax)


    # G + LB
    center_G_LB = [[0, 12.5/10], [1/10, 12.5/10], [2/10, 12.5/10], [2.5/10, 12.5/10], [3/10, 12.5/10]]
    edge_G_LB = [
        [[0, 4.1531], [0, 22.1902]],
        [[1.3039, 4.2378], [0.6576, 21.9206]],
        [[2.6128, 4.3159], [1.3177, 21.7441]],
        [[3.2733, 4.3624], [1.6440, 21.6433]],
        [[3.9437, 4.4311], [1.9665, 21.5689]]
    ]
    plot_tie_lines(ax, center_G_LB, edge_G_LB, c="tab:grey")
    fill_color(edge_G_LB, ax)
    
    # LB + LAB
    center_LB_LAB = [[10, 45], [10, 40], [10, 38]]
    edge_LB_LAB = [
        [[4.3168, 31.2870], [21.8533, 73.9867]],
        [[4.4412, 28.2210], [21.3746, 64.7041]],
        [[4.5254, 26.9304], [21.2044, 61.2696]]
    ]
    plot_tie_lines(ax, center_LB_LAB, edge_LB_LAB, c="tab:grey")
    fill_color([[[3, 75], [22, 76]]] + edge_LB_LAB, ax, c="tab:blue")

    # LAB + LA
    center_LAB_LA = [
        # [130/4, 40/4],
        [120/4, 50/4],
        [130/4, 50/4],
        [140/4, 50/4],
        [160/4, 55/4],
        # [45, 15],
        [47, 19],
        [60, 20]]
    edge_LAB_LA = [
        # [[90.3333/4, 70.2929/4], [152.7865/4, 22.8993/4]],
        [[91.1910/4, 71.7885/4], [156.7620/4, 22.5345/4]], 
        [[97.4319/4, 72.6223/4], [171.5426/4, 21.3663/4]],
        [[103.2168/4, 73.9762/4], [185.9900/4, 20.3843/4]],
        [[118.4965/4, 78.5602/4], [225.7040/4, 18.3693/4]],
        # [[36.65, 19.], [66.4, 4.578]],
        [[41.6492, 21.3472], [88.2805, 4.2819]],
        [[53.6089, 22.2076], [108.7682, 4.1077]]
    ]
    plot_tie_lines(ax, center_LAB_LA, edge_LAB_LA, c="tab:grey")
    fill_color(edge_LAB_LA + [[[108, 24], [108.7682, 4.1077]]], ax, c="tab:blue")

    # DC + LA
    center_LA_DC = [
        [16.6, 13.6],
        [18, 14.5],
        [19.5, 15.5],
        [21, 16.5],
        [52/2, 28/2]
        # [104/4, 56/4]
    ]
    edge_LA_DC = [
        [[21.6557, 2.9274], [13.2678, 21.2943]],
        [[24.7437, 3.4352], [14.2187, 21.0731]],
        [[27.2349, 3.8950], [15.9725, 21.1980]],
        [[30.2194, 4.6432], [17.8965, 20.8375]],
        [[67.2622/2, 10.1195/2], [40.1760/2, 41.6459/2]]
        # [[134.0522/4, 19.6533/4], [79.7745/4, 87.2652/4]]
    ]
    plot_tie_lines(ax, center_LA_DC, edge_LA_DC, "s")
    edge = [[edge_3p[2], edge_3p[0]]] + edge_LA_DC + [edge_LAB_LA[0][::-1]]
    fill_color(edge, ax, c="tab:orange")

    # DC + G + LB
    edge_G_LB_DC = [[3.9437, 4.4311], [1.9665, 21.5689], [7.1505, 21.9526]]
    x = [i[0] for i in edge_G_LB_DC]
    y = [i[1] for i in edge_G_LB_DC]
    ax.fill(x, y, c="tab:purple", alpha=0.5)

    # LB + DC
    edge_LB_DC = [
        [edge_G_LB_DC[1], edge_G_LB_DC[2]],
        [[3.7, 22.5], [9, 22.5]],
        [[4.1, 23.0], [10, 24.5]],
        [[4.65, 24.0], [14.5, 34]],
        [[4.6, 25.0], [16.5, 40]],
        [[4.55, 26.0], [19.0, 50]],
        edge_LB_LAB[-1]
    ]
    fill_color(edge_LB_DC, ax, c="tab:brown")

    # DC
    edge_DC = [p[1] for p in edge_LB_DC[::-1]] + [p[1] for p in edge_G_DC[::-1]] + [p[1] for p in edge_LA_DC]
    edge_DC.append(edge_LAB_LA[0][0])
    edge_DC += [
        [23.5, 20],
        [22, 22.5],
        [21, 25],
        [20.5, 30],
        edge_LB_LAB[-1][1]
    ]
    x = [i[0] for i in edge_DC]
    y = [i[1] for i in edge_DC]
    ax.fill(x, y, c="tab:green", alpha=0.5)

    # ax.plot(10, 24, 'H', c="tab:green")
    # ax.plot(10, 23, 'H', c="tab:green")
    # ax.plot(12.5, 25, 'H', c="tab:green")
    # ax.plot(15, 25, 'H', c="tab:green")
    # ax.plot(10, 25, 'p', c="tab:brown")
    # ax.plot(10, 26, 'p', c="tab:brown")
    # ax.plot(10, 28, 'p', c="tab:brown")
    # ax.plot(10, 34, "o", c="tab:red")
    # ax.plot(4.5, 14.7, "<", fillstyle="none", c="tab:green")
    # ax.plot(3, 6, "x", c="tab:green")
    # ax.plot(5.5, 3.85, "x", c="tab:red")
    # ax.plot(6, 3, "x", c="tab:blue")
    # ax.plot(20, 20, "x", c="tab:pink")
    # ax.plot(2.5, 22.5, "o", c="tab:grey", fillstyle="none")
    # ax.plot(3.5, 22.5, "o", c="tab:grey", fillstyle="none")
    # ax.plot(4.0, 22.5, "p", c="tab:brown")
    # ax.plot(4.5, 22.5, "p", c="tab:brown")
    # ax.plot(5, 22.5, "p", c="tab:brown")
    # ax.plot(3.5, 23, "o", c="tab:grey", fillstyle="none")
    # ax.plot(4, 23, "o", c="tab:grey", fillstyle="none")
    # ax.plot(4.5, 23, "p", c="tab:brown")
    # ax.plot(4, 24, "o", c="tab:grey", fillstyle="none")
    # ax.plot(4.5, 24, "o", c="tab:grey", fillstyle="none")
    # ax.plot(5, 24, "p", c="tab:brown")
    # ax.plot(2.5, 25, "o", c="tab:grey", fillstyle="none")
    # ax.plot(4.5, 25, "o", c="tab:grey", fillstyle="none")
    # ax.plot(5, 25, "p", c="tab:brown")
    # ax.plot(4.5, 26, "o", c="tab:grey", fillstyle="none")
    # ax.plot(5, 26, "p", c="tab:brown")
    # ax.plot(25, 20, "o", c="tab:grey", fillstyle="none")
    # ax.plot(22.5, 20, "H", c="tab:green")
    # ax.plot(22.5, 22.5, "o", c="tab:grey", fillstyle="none")
    # ax.plot(20, 22.5, "H", c="tab:green")
    # ax.plot(22.5, 25, "o", c="tab:grey", fillstyle="none")
    # ax.plot(20, 25, "H", c="tab:green")
    # ax.plot(22.5, 30, "o", c="tab:grey", fillstyle="none")
    # ax.plot(20, 30, "H", c="tab:green")
    # ax.plot(15, 34, "H", c="tab:green")
    # ax.plot(12.5, 34, "p", c="tab:brown")
    # ax.plot(11, 34, "o", c="tab:red")
    # ax.plot(17.5, 40, "H", c="tab:green")
    # ax.plot(15, 40, "p", c="tab:brown")
    # ax.plot(12.5, 40, "o", c="tab:red")
    # # plt.plot(20, 60, "o", c="tab:blue")
    # ax.plot(22.5, 60, "o", c="tab:grey", fillstyle="none")
    # # plt.plot(22.5, 50, "o", c="tab:grey", fillstyle="none")

    # ax.plot(6, 30, "x")
    # ax.plot(5, 30, "x")

    # ax.plot(15, 50, "x")
    # ax.plot(15, 40, "x")

    # ax.plot(110/4, 50/4, "o")
    # ax.plot(115/4, 50/4, "*")
    # ax.plot(116/4, 50/4, "*")

    # plt.plot(130/4, 50/4, "p")
    # ax.plot(2.5, 80/4, "p")

    # ax.set_xlim(0, 110)
    # ax.set_ylim(0, 75)

    ax.set_xlim(0, 70)
    ax.set_ylim(0, 70)
    ax.set_xlabel(r"$\phi_A$", fontsize="x-large")
    ax.set_ylabel(r"$\phi_B$", fontsize='x-large')
    # ax.set_title(r"$D_r=0.1, \eta=-2,\eta_{AB}=-\eta_{BA}=0.5, \rho_0=\bar{\rho}_{A,B}=10$", fontsize="x-large")

    # ax.text(40, 40, "LAB", fontsize="x-large")
    # ax.text(15, 30, "DC", fontsize="x-large")
    # ax.text(40, 2, "LA", fontsize="x-large")
    # ax.text(0.5, 0.5, "G", fontsize="x-large")
    # ax.text(0.5, 40, "LB", fontsize="x-large")
    # ax.text(10, 60, "LB+LAB", fontsize="x-large")
    # ax.text(0.5, 8, "LB+G", fontsize="x-large", rotation=-90)
    # ax.text(8, 30, "LB+DC", fontsize="x-large")
    # ax.text(3, 15, "LB+G+DC", fontsize="x-large", rotation=-90)
    # ax.text(7.5, 18, "DC+G", fontsize="x-large", rotation=0)
    # ax.text(7, 4, "LA+G+DC", fontsize="x-large", rotation=0)
    # ax.text(14, 0.5, "LA+G", fontsize="x-large", rotation=0)
    # ax.text(24, 10, "LA+DC", fontsize="x-large", rotation=0)
    # ax.text(50, 15, "LAB+LA", fontsize="x-large", rotation=0)




if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1, figsize=(12, 12), constrained_layout=True)
    binodal(ax)
    plt.show()
    plt.close()
    
        
