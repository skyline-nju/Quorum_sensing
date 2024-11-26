import numpy as np
import matplotlib.pyplot as plt


root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"

def plot_colored_line(x, y, c, ax=None, cmap='viridis'):
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm
    if ax is None:
        ax = plt.gca()

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(c.min(), c.max())
    # norm = plt.Normalize(0, 1)

    lc = LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping
    lc.set_array(c)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    # plt.colorbar(line, ax=ax)

def plot_colourline(x,y,c, ax=None):
    import matplotlib.cm as cm
    c = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
    if ax is None:
        ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i])
    return


if __name__ == "__main__":
    dx = 2
    prefix = f"{root_rudabeh}/QS_MPI_noself/L80_r20_20_20/coarse_grain_dx{dx:d}"
    basename = "L80_80_Dr0.100_k0.70_p40_40_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"

    fname = f"{prefix}/{basename}"

    dt = 0.01
    with np.load(fname, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    
    # beg = 19
    beg = 119
    rho_A = fields[beg:, 0, 0, 0]
    rho_B = fields[beg:, 1, 0, 0]
    theta_A = np.arctan2(fields[beg:, 4, 0, 0], fields[beg:, 2, 0, 0])
    theta_B = np.arctan2(fields[beg:, 5, 0, 0], fields[beg:, 3, 0, 0])

    t = t[beg:]

    fig = plt.figure(constrained_layout=True, figsize=(12, 4.2))
    subfigs = fig.subfigures(1, 2, width_ratios=[2, 1])
    # subfigs = fig.subfigures(2, 1, height_ratios=[0.35, 0.65])

    (ax1, ax2) = subfigs[0].subplots(2, 1, sharex=True)
    ax1.plot(t, rho_A)
    ax1.plot(t, rho_B)
    ax2.plot(t, theta_A, "-", ms=4, fillstyle="none")
    ax2.plot(t, theta_B, "-", ms=4, fillstyle="none")
    ax1.set_xlim(t[0], t[-1])
    ax3 = subfigs[1].subplots(1, 1)
    plot_colourline(rho_A, rho_B, t, ax=ax3)
    # plot_colored_line(rho_A, rho_B, t, ax=ax3)
    ax1.set_ylabel(r"$\rho_{A,B}(\mathbf{r},t)$", fontsize="large")
    ax2.set_ylabel(r"$\theta_{A,B}(\mathbf{r},t)$", fontsize="large")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax3.set_xlabel(r"$\rho_A$", fontsize="large")
    ax3.set_ylabel(r"$\rho_B$", fontsize="large")
    ax1.set_title("(a)")
    ax2.set_title("(b)")
    ax3.set_title("(c)")
    
    title = r"$D_r=0.1, \eta=-2,\eta_{AB}=-\eta_{BA}=0.5, \phi_A=\phi_B=40,\rho_0=\bar{\rho}_{A,B}=20$"
    fig.suptitle(title)
    plt.show()
    plt.close()
