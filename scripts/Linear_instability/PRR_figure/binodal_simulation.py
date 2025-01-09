import numpy as np
import matplotlib.pyplot as plt

# from profile import plot_time_ave_density_profile_SP


def get_slope(p1, p2):
    return (p1[1] - p2[1]) / (p1[0] - p2[0])

def plot_tri(ax, edge, center, c):
    for i in range(len(edge)):
        # ax.plot([edge[i-1][0], edge[i][0]], [edge[i-1][1], edge[i][1]], c=c, linestyle=linestyle)
        ax.plot([edge[i][0], center[0]], [edge[i][1], center[1]], c=c, linestyle="dotted")

def plot_bi(ax, edge, c, linestyle="dotted"):
    ax.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], c=c, linestyle=linestyle)


def plot_tie_lines(ax, centers, edges, marker="o", c=None, ms=None):
    for i, center in enumerate(centers):
        x = [j[0] for j in edges[i]]
        y = [j[1] for j in edges[i]]
        # line, = ax.plot(center[0], center[1], "s", c=c, fillstyle="none")
        ax.plot(x, y, marker, fillstyle="none", c=c, ms=ms)
        plot_bi(ax, edges[i], c=c)


def fill_color(edge, ax, c="tab:blue", alpha=0.25):
    n = len(edge)
    x, y = [], []
    for p in edge:
        x.append(p[0][0])
        y.append(p[0][1])
    for p in edge[::-1]:
        x.append(p[1][0])
        y.append(p[1][1])
    ax.fill(x, y, c=c, alpha=alpha)


def plot_PD(ax=None, show_tie_line=True, label_font_size="xx-large", ms=3):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    
    # LB + DC + G
    center_3p = np.array([15, 10]) / 10
    edge_3p = np.array([[13.215, 21.3376], [4.4581, 1.1997], [20.6741, 2.8125]]) / 10
    x = [i[0] for i in edge_3p]
    y = [i[1] for i in edge_3p]
    # line, = ax.plot(center_3p[0], center_3p[1], ">")
    # ax.plot(x, y, ">", fillstyle="none", c=line.get_c())
    if show_tie_line:
        plot_tri(ax, edge_3p, center_3p, c="tab:grey")

    x = [i[0] for i in edge_3p]
    y = [i[1] for i in edge_3p]
    ax.fill(x, y, c="tab:orange", alpha=0.25)

    # DC + G
    center_G_DC = np.array([[9.5, 13.075], [8.625, 13.36], [7.5, 13.725], [6.5, 14.05]]) / 10
    edge_G_DC = np.array([[[4.1251, 1.2010], [11.8641, 21.5710]],
                 [[3.9359, 1.1988], [10.4761, 21.9334]],
                 [[3.6668, 1.2039], [8.644, 22.1592]],
                 [[3.404, 1.3104], [7.1505, 21.9526]]
                ])/10
    if show_tie_line:
        plot_tie_lines(ax, center_G_DC, edge_G_DC, "o", c="tab:grey", ms=ms)
    edge = np.vstack((np.array([[edge_3p[1], edge_3p[0]]]), edge_G_DC))
    fill_color(edge, ax, c="tab:cyan")

    # G + LA
    center_G_LA = np.array([
        [12.5, 0],
        [12.5, 0.5],
        [12.5, 1],
        [12.5, 2],
        # [12.5, 3]
    ])/10
    edge_G_LA = np.array([
        [[4.1229, 0], [22.1689, 0]],
        [[4.1984, 0.2967], [21.8290, 0.7313]],
        [[4.2891, 0.5983], [21.5045, 1.4470]],
        [[4.4534, 1.2062], [20.6741, 2.8125]],
        # [[4.6308, 1.7986], [19.6164, 4.1027]]
    ]) / 10
    if show_tie_line:
        plot_tie_lines(ax, center_G_LA, edge_G_LA, c="tab:grey", ms=ms)
    fill_color(edge_G_LA, ax)


    # G + LB
    center_G_LB = np.array([[0, 12.5], [1, 12.5], [2, 12.5], [2.5, 12.5], [3, 12.5]]) / 10
    edge_G_LB = np.array([
        [[0, 4.1531], [0, 22.1902]],
        [[1.3039, 4.2378], [0.6576, 21.9206]],
        [[2.6128, 4.3159], [1.3177, 21.7441]],
        [[3.2733, 4.3624], [1.6440, 21.6433]],
        [[3.9437, 4.4311], [1.9665, 21.5689]]
    ]) / 10
    if show_tie_line:
        plot_tie_lines(ax, center_G_LB, edge_G_LB, c="tab:grey", ms=ms)
    fill_color(edge_G_LB, ax)
    
    # LB + LAB
    center_LB_LAB = np.array([[10, 45], [10, 40], [10, 38]]) / 10
    edge_LB_LAB = np.array([
        [[4.3168, 31.2870], [21.8533, 73.9867]],
        [[4.4412, 28.2210], [21.3746, 64.7041]],
        [[4.5254, 26.9304], [21.2044, 61.2696]]
    ]) / 10
    if show_tie_line:
        plot_tie_lines(ax, center_LB_LAB, edge_LB_LAB, c="tab:grey", ms=ms)
    edge = np.vstack((np.array([[[3, 75], [22, 76]]])/10, edge_LB_LAB))
    fill_color(edge, ax, c="tab:blue")

    # LAB + LA
    center_LAB_LA = np.array([
        # [130/4, 40/4],
        [120/4, 50/4],
        [130/4, 50/4],
        [140/4, 50/4],
        [160/4, 55/4],
        # [45, 15],
        [47, 19],
        [60, 20]]) / 10
    edge_LAB_LA = np.array([
        # [[90.3333/4, 70.2929/4], [152.7865/4, 22.8993/4]],
        [[91.1910/4, 71.7885/4], [156.7620/4, 22.5345/4]], 
        [[97.4319/4, 72.6223/4], [171.5426/4, 21.3663/4]],
        [[103.2168/4, 73.9762/4], [185.9900/4, 20.3843/4]],
        [[118.4965/4, 78.5602/4], [225.7040/4, 18.3693/4]],
        # [[36.65, 19.], [66.4, 4.578]],
        [[41.6492, 21.3472], [88.2805, 4.2819]],
        [[53.6089, 22.2076], [108.7682, 4.1077]]
    ]) / 10
    if show_tie_line:
        plot_tie_lines(ax, center_LAB_LA, edge_LAB_LA, c="tab:grey", ms=ms)
    edge = np.vstack((edge_LAB_LA, np.array([[[108, 24], [108.7682, 4.1077]]])/10))
    fill_color(edge, ax, c="tab:blue")

    # DC + LA
    center_LA_DC = np.array([
        [16.6, 13.6],
        [18, 14.5],
        [19.5, 15.5],
        [21, 16.5],
        [52/2, 28/2]
        # [104/4, 56/4]
    ]) / 10
    edge_LA_DC = np.array([
        [[21.6557, 2.9274], [13.2678, 21.2943]],
        [[24.7437, 3.4352], [14.2187, 21.0731]],
        [[27.2349, 3.8950], [15.9725, 21.1980]],
        [[30.2194, 4.6432], [17.8965, 20.8375]],
        [[67.2622/2, 10.1195/2], [40.1760/2, 41.6459/2]]
        # [[134.0522/4, 19.6533/4], [79.7745/4, 87.2652/4]]
    ]) / 10
    if show_tie_line:
        plot_tie_lines(ax, center_LA_DC, edge_LA_DC, c="tab:grey", ms=ms)
    edge = np.vstack((np.array([[edge_3p[2], edge_3p[0]]]), edge_LA_DC))
    edge = np.vstack((edge, np.array([edge_LAB_LA[0][::-1]])))
    fill_color(edge, ax, c="tab:cyan")

    # DC + G + LB
    edge_G_LB_DC = np.array([[3.9437, 4.4311], [1.9665, 21.5689], [7.1505, 21.9526]]) / 10
    x = [i[0] for i in edge_G_LB_DC]
    y = [i[1] for i in edge_G_LB_DC]
    ax.fill(x, y, c="tab:orange", alpha=0.25)

    # LB + DC
    edge_LB_DC = np.array([
        [edge_G_LB_DC[1], edge_G_LB_DC[2]],
        [[3.7/10, 22.5/10], [9/10, 22.5/10]],
        [[4.1/10, 23.0/10], [1, 2.45]],
        [[4.65/10, 2.4], [1.45, 3.4]],
        [[4.6/10, 2.5], [1.65, 4]],
        [[4.55/10, 2.6], [1.9, 5]],
        edge_LB_LAB[-1]
    ]) 
    fill_color(edge_LB_DC, ax, c="tab:cyan")

    # DC
    edge_DC = [p[1] for p in edge_LB_DC[::-1]] + [p[1] for p in edge_G_DC[::-1]] + [p[1] for p in edge_LA_DC]
    edge_DC.append(edge_LAB_LA[0][0])
    edge_DC += [
        [22.5/10, 20/10],
        [22.4/10, 22.5/10],
        [22.3/10, 25/10],
        [22.2/10, 30/10],
        edge_LB_LAB[-1][1]
    ]
    x = [i[0] for i in edge_DC]
    y = [i[1] for i in edge_DC]
    ax.fill(x, y, c="tab:green", alpha=0.3)

    # plt.plot(10, 23, 'H', c="tab:green")
    # plt.plot(10, 24, 'H', c="tab:green")
    # plt.plot(12.5, 25, 'H', c="tab:green")
    # plt.plot(15, 25, 'H', c="tab:green")

    # plt.plot(10, 25, 'p', c="tab:brown")
    # plt.plot(10, 26, 'p', c="tab:brown")
    # plt.plot(10, 28, 'p', c="tab:brown")

    # plt.plot(10, 34, "o", c="tab:red")

    # plt.plot(4.5, 14.7, "<", fillstyle="none", c="tab:green")
    # # plt.plot(4, 12.5, "<", fillstyle="none", c="tab:green")


    # plt.plot(3, 6, "x", c="tab:green")
    # plt.plot(5.5, 3.85, "x", c="tab:red")
    # plt.plot(6, 3, "x", c="tab:blue")
    # # plt.plot(6.5, 1.6, "x", c="tab:blue")

    # plt.plot(20, 20, "x", c="tab:pink")

    # plt.plot(2.5, 22.5, "o", c="tab:grey", fillstyle="none")
    # plt.plot(3.5, 22.5, "o", c="tab:grey", fillstyle="none")
    # plt.plot(4.0, 22.5, "p", c="tab:brown")
    # plt.plot(4.5, 22.5, "p", c="tab:brown")
    # plt.plot(5, 22.5, "p", c="tab:brown")

    # plt.plot(3.5, 23, "o", c="tab:grey", fillstyle="none")
    # plt.plot(4, 23, "o", c="tab:grey", fillstyle="none")
    # plt.plot(4.5, 23, "p", c="tab:brown")

    # plt.plot(4, 24, "o", c="tab:grey", fillstyle="none")
    # plt.plot(4.5, 24, "o", c="tab:grey", fillstyle="none")
    # plt.plot(5, 24, "p", c="tab:brown")

    # plt.plot(2.5, 25, "o", c="tab:grey", fillstyle="none")
    # plt.plot(4.5, 25, "o", c="tab:grey", fillstyle="none")
    # plt.plot(5, 25, "p", c="tab:brown")
    # plt.plot(4.5, 26, "o", c="tab:grey", fillstyle="none")
    # plt.plot(5, 26, "p", c="tab:brown")
    # plt.plot(6, 30, "p", c="tab:brown")


    # plt.plot(25, 20, "o", c="tab:grey", fillstyle="none")
    # plt.plot(22.5, 20, "H", c="tab:green")
    # plt.plot(22.5, 22.5, "o", c="tab:grey", fillstyle="none")
    # plt.plot(20, 22.5, "H", c="tab:green")
    # plt.plot(22.5, 25, "o", c="tab:grey", fillstyle="none")
    # plt.plot(20, 25, "H", c="tab:green")
    # plt.plot(22.5, 30, "o", c="tab:grey", fillstyle="none")
    # plt.plot(20, 30, "H", c="tab:green")

    # plt.plot(22.5, 34, "o", c="tab:grey", fillstyle="none")
    # plt.plot(20, 34, "H", c="tab:green")



    # plt.plot(15, 34, "H", c="tab:green")
    # plt.plot(12.5, 34, "p", c="tab:brown")
    # plt.plot(11, 34, "o", c="tab:red")
    # plt.plot(17.5, 40, "H", c="tab:green")
    # plt.plot(15, 40, "p", c="tab:brown")
    # plt.plot(12.5, 40, "o", c="tab:red")
    # # plt.plot(20, 60, "o", c="tab:blue")
    # plt.plot(22.5, 60, "o", c="tab:grey", fillstyle="none")
    # # plt.plot(22.5, 50, "o", c="tab:grey", fillstyle="none")

    # # plt.plot(4.4, 23.6, "*", c="k")
    # # plt.plot(21.8, 56.5, "*", c="k")
    # # plt.plot(10, 34, "*", c="k")
    # # plt.plot(12.5, 10, "x", c="k")

    # # plt.plot(3, 17.3, "o")
    # plt.plot(4, 7.9, "o")



    # # plt.plot(6, 30, "x")
    # # plt.plot(5, 30, "x")

    # # plt.plot(15, 50, "x")
    # # plt.plot(15, 40, "x")

    # # plt.plot(110/4, 50/4, "o")
    # # plt.plot(115/4, 50/4, "*")
    # # plt.plot(116/4, 50/4, "*")

    # # plt.plot(130/4, 50/4, "p")
    # # plt.plot(4, 8, "p")
    # # plt.plot(4.3, 8, "p")


    # # ax.set_xlim(0, 110)
    # # ax.set_ylim(0, 75)

    if show_tie_line:
        ax.set_xlim(0, 4.4)
        ax.set_ylim(0, 6.5)
    # ax.set_xlabel(r"$\bar{\rho}_A/\rho_0$", fontsize="large")
    # ax.set_ylabel(r"$\bar{\rho}_B/\rho_0$", fontsize='large')
    # ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
    # ax.set_xticklabels([0, 1, 2, 3, 4, 5, 6, 7])

    # ax.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
    # ax.set_yticklabels([0, 1, 2, 3, 4, 5, 6, 7])
    # ax.set_title(r"$D_r=0.1, \eta=-2,\eta_{AB}=-\eta_{BA}=0.5$", fontsize="x-large")

    if show_tie_line:
        fs ="x-large"
        ax.text(3.5, 3, "LAB", fontsize=fs)
        ax.text(1.5, 3, "CCB", fontsize=fs)
        ax.text(3, 0.2, "LA", fontsize=fs)
        ax.text(0.05, 0.05, "G", fontsize=fs)
        ax.text(0.05, 3, "LB", fontsize=fs)

        ax.text(0.5, 4, "LB+LAB", fontsize=fs, rotation=70)
        ax.text(0.05, 0.8, "LB+G", fontsize=fs, rotation=-90)
        ax.text(1, 3, "LB+CCB", fontsize=fs, rotation=60)
        ax.text(0.3, 1.1, "LB+G+CCB", fontsize=fs, rotation=-90)
        ax.text(0.7, 1.4, "G+CCB", fontsize=fs, rotation=70)
        ax.text(0.7, 0.4, "G+LA+CCB", fontsize=fs, rotation=0)
        ax.text(1.4, 0.05, "LA+G", fontsize=fs, rotation=0)
        ax.text(2, 1.0, "CCB+LA", fontsize=fs, rotation=-50)
        ax.text(3, 1.4, "LAB+LA", fontsize=fs, rotation=-30)

        ax.text(3.8, 0.1, r"$\bar{\rho}_A/\rho_0$", fontsize=label_font_size)
        ax.text(0.05, 5.5, r"$\bar{\rho}_B/\rho_0$", fontsize=label_font_size, rotation=90)

    # plot_PD_composition(state, extent, ax=ax, fill=False, set_xy_lims=False, scale_factor=10)

    if flag_show:
        plt.show()
        plt.close()
    

if __name__ == "__main__":
    plot_PD()
