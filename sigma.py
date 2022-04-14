import numpy as np
import matplotlib.pyplot as plt
from phase_diagram import read_txt

if __name__ == "__main__":
    Lx=20
    Ly=5
    phiA=80
    phiB=80
    rho0=80
    k=0.7
    eta=0


    fout = f"fig/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_e{eta:g}.png"
    fnpz = f"data/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_e{eta:g}.npz"
    fdat = f"data/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_e{eta:g}.dat"
    paras = read_txt(fdat, ["Dr", "alpha", "state"])
    Dr = paras["Dr"]
    alpha = paras["alpha"]
    state = paras["state"]

    with np.load(fnpz, "r") as data:
        # alpha = data["alpha"]
        # Dr = data["Dr"]
        v_std = data["v_std"]
    

    for D in [0.01, 0.1, 1, 3]:
        mask = Dr == D
        plt.plot(alpha[mask], v_std[mask], "-o", label=r"$D_r=%g$" % D)
    
    # plt.xscale("log")
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\sigma_v$")
    plt.legend()
    plt.show()
    plt.close()

