from types import DynamicClassAttribute
import numpy as np
import matplotlib.pyplot as plt
from corr_length import get_Lx, get_para


root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def ori_ave(C_xy, r_radial, r_module):
    C_r = np.zeros(r_radial.size)
    half_dr = (r_radial[1] - r_radial[0]) / 2

    ny, nx = C_xy.shape
    for i, r in enumerate(r_radial):
        if i == 0:
            C_r[0] = C_xy[ny//2, nx//2]
        else:
            mask = np.logical_and(r_module >= r - half_dr, r_module < r + half_dr)
            C_r[i] = np.mean(C_xy[mask])
    return C_r


def get_suptitle(fname, t):
    para = get_para(fname)
    title = r"$L=%g,D_r=%g,\eta=%g,\eta_{AB}=-\eta_{BA}=%g,\rho_0=\bar{\rho}_{A,B}=%g,\phi_A=%g,\phi_B=%g,t=%g$" % (
        para['Lx'], para["Dr"], para["eta"], para["eta_AB"], para["rho0"], para['pA'], para['pB'], t * para['h']
    )
    return title


def plot_corr_func(fin, dx=1, beg=None, end=None, title=None):
    with np.load(fin, "r") as data:
        fields = data["fields"]
        x, y, t = data["x"], data["y"], data["t"]
    if beg is None:
        beg = 0
    if end is None:
        end = t.size
    r_module = None
    r_radial = None
    for i in range(beg, end):
        print("cal corr func for %d/%d frame, t=%g" % (i, t.size, t[i]))
        rhoA_m, rhoB_m = np.mean(fields[i][0]), np.mean(fields[i][1])
        rhoA_q = np.fft.fft2(fields[i][0]-rhoA_m, norm="ortho")
        rhoB_q = np.fft.fft2(fields[i][1]-rhoB_m, norm="ortho")
        S_rho_AA = rhoA_q * rhoA_q.conj()
        S_rho_BB = rhoB_q * rhoB_q.conj()
        S_rho_AB = rhoA_q * rhoB_q.conj()
        C_rho_AA = np.fft.ifft2(S_rho_AA).real # C_rho_AA[0, 0] = <(rho_A - rho_A_m)^2>
        C_rho_BB = np.fft.ifft2(S_rho_BB).real
        C_rho_AB = np.fft.ifft2(S_rho_AB).real
        C_rho_AA = np.fft.fftshift(C_rho_AA)
        C_rho_BB = np.fft.fftshift(C_rho_BB)
        C_rho_AB = np.fft.fftshift(C_rho_AB)

        mx_A, mx_B, my_A, my_B = fields[i][2:6]
        print("rho:", rhoA_m, rhoB_m, "m:", np.mean(mx_A), np.mean(my_A), np.mean(mx_B), np.mean(my_B))
        mx_A_q = np.fft.fft2(mx_A, norm="ortho")
        mx_B_q = np.fft.fft2(mx_B, norm="ortho")
        my_A_q = np.fft.fft2(my_A, norm="ortho")
        my_B_q = np.fft.fft2(my_B, norm="ortho")
        C_m_AA = np.fft.fftshift(np.fft.ifft2(mx_A_q * mx_A_q.conj() + my_A_q * my_A_q.conj()).real)
        C_m_BB = np.fft.fftshift(np.fft.ifft2(mx_B_q * mx_B_q.conj() + my_B_q * my_B_q.conj()).real)
        C_m_AB = np.fft.fftshift(np.fft.ifft2(mx_A_q * mx_B_q.conj() + my_A_q * my_B_q.conj()).real)

        mask_A = fields[i][0] > 0
        mask_B = fields[i][1] > 0
        px_A, py_A, px_B, py_B = np.zeros((4, mx_A.shape[0], mx_A.shape[1]))
        px_A[mask_A] = mx_A[mask_A] / fields[i][0][mask_A]
        py_A[mask_A] = my_A[mask_A] / fields[i][0][mask_A]
        px_B[mask_B] = mx_B[mask_B] / fields[i][1][mask_B]
        py_B[mask_B] = my_B[mask_B] / fields[i][1][mask_B]
        px_A_q = np.fft.fft2(px_A, norm="ortho")
        px_B_q = np.fft.fft2(px_B, norm="ortho")
        py_A_q = np.fft.fft2(py_A, norm="ortho")
        py_B_q = np.fft.fft2(py_B, norm="ortho")

        S_p_AA = px_A_q * px_A_q.conj() + py_A_q * py_A_q.conj()
        S_p_BB = px_B_q * px_B_q.conj() + py_B_q * py_B_q.conj()
        S_p_AB = px_A_q * px_B_q.conj() + py_A_q * py_B_q.conj()
        C_p_AA = np.fft.fftshift(np.fft.ifft2(S_p_AA).real)
        C_p_BB = np.fft.fftshift(np.fft.ifft2(S_p_BB).real)
        C_p_AB = np.fft.fftshift(np.fft.ifft2(S_p_AB).real)

        vx_A, vx_B, vy_A, vy_B = fields[i][8:12]
        vx_A_q = np.fft.fft2(vx_A, norm="ortho")
        vx_B_q = np.fft.fft2(vx_B, norm="ortho")
        vy_A_q = np.fft.fft2(vy_A, norm="ortho")
        vy_B_q = np.fft.fft2(vy_B, norm="ortho")

        S_v_AA = vx_A_q * vx_A_q.conj() + vy_A_q * vy_A_q.conj()
        S_v_BB = vx_B_q * vx_B_q.conj() + vy_B_q * vy_B_q.conj()
        S_v_AB = vx_A_q * vx_B_q.conj() + vy_A_q * vy_B_q.conj()
        C_v_AA = np.fft.fftshift(np.fft.ifft2(S_v_AA).real)
        C_v_BB = np.fft.fftshift(np.fft.ifft2(S_v_BB).real)
        C_v_AB = np.fft.fftshift(np.fft.ifft2(S_v_AB).real)

        s_A, s_B = fields[i][6:8]
        s_A_q = np.fft.fft2(s_A-np.mean(s_A), norm="ortho")
        s_B_q = np.fft.fft2(s_B-np.mean(s_B), norm="ortho")

        S_s_AA = s_A_q * s_A_q.conj()
        S_s_BB = s_B_q * s_B_q.conj()
        S_s_AB = s_A_q * s_B_q.conj()
        C_s_AA = np.fft.fftshift(np.fft.ifft2(S_s_AA).real)
        C_s_BB = np.fft.fftshift(np.fft.ifft2(S_s_BB).real)
        C_s_AB = np.fft.fftshift(np.fft.ifft2(S_s_AB).real)


        # print(np.min(fields[i][6:8]), np.max(fields[i][6:8]))
        print(np.mean(vx_A * vx_B + vy_A * vy_B))
        if r_module is None:
            ny, nx = C_rho_AA.shape
            X = (np.arange(nx) - nx // 2) * dx
            Y = (np.arange(ny) - ny // 2) * dx
            xx, yy = np.meshgrid(X, Y)
            r_module = np.sqrt(xx ** 2 + yy ** 2)
            r_radial = np.arange(nx // 2) * dx

        Cr_rho_AA = ori_ave(C_rho_AA, r_radial, r_module)
        Cr_rho_BB = ori_ave(C_rho_BB, r_radial, r_module)
        Cr_rho_AB = ori_ave(C_rho_AB, r_radial, r_module)
        Cr_m_AA = ori_ave(C_m_AA, r_radial, r_module)
        Cr_m_BB = ori_ave(C_m_BB, r_radial, r_module)
        Cr_m_AB = ori_ave(C_m_AB, r_radial, r_module)
        Cr_p_AA = ori_ave(C_p_AA, r_radial, r_module)
        Cr_p_BB = ori_ave(C_p_BB, r_radial, r_module)
        Cr_p_AB = ori_ave(C_p_AB, r_radial, r_module)
        Cr_v_AA = ori_ave(C_v_AA, r_radial, r_module)
        Cr_v_BB = ori_ave(C_v_BB, r_radial, r_module)
        Cr_v_AB = ori_ave(C_v_AB, r_radial, r_module)
        Cr_s_AA = ori_ave(C_s_AA, r_radial, r_module)
        Cr_s_BB = ori_ave(C_s_BB, r_radial, r_module)
        Cr_s_AB = ori_ave(C_s_AB, r_radial, r_module)
        fig, axes = plt.subplots(2, 2, figsize=(8, 6.2), constrained_layout=True, sharex=True)

        for ax in axes.flat:
            ax.axhline(0, linestyle="dotted", c="tab:grey")

        axes[0, 0].plot(r_radial, Cr_rho_AA, label="AA")
        axes[0, 0].plot(r_radial, Cr_rho_BB, label="BB")
        axes[0, 0].plot(r_radial, Cr_rho_AB, label="AB")
        axes[0, 1].plot(r_radial, Cr_p_AA, label="AA")
        axes[0, 1].plot(r_radial, Cr_p_BB, label="BB")
        axes[0, 1].plot(r_radial, Cr_p_AB, label="AB")
        axes[1, 0].plot(r_radial, Cr_s_AA, label="AA")
        axes[1, 0].plot(r_radial, Cr_s_BB, label="BB")
        axes[1, 0].plot(r_radial, Cr_s_AB, label="AB")
        axes[1, 1].plot(r_radial, Cr_v_AA, label="AA")
        axes[1, 1].plot(r_radial, Cr_v_BB, label="BB")
        axes[1, 1].plot(r_radial, Cr_v_AB, label="AB")
        # ax1.set_xlim(0, L/3)

        titles = [r"(a) $\langle \delta \rho_S(0) \delta \rho_{S'}(\mathbf{r})\rangle_r$",
                    r"(b) $\langle \mathbf{p}_S(0)\cdot\mathbf{p}_{S'}(\mathbf{r})\rangle_r$",
                    r"(c) $\langle \delta\bar{v}_S(0)\delta\bar{v}_{S'}(\mathbf{r})\rangle_r$",
                    r"(d) $\langle \mathbf{v}_S(0)\cdot\mathbf{v}_{S'}(\mathbf{r})\rangle_r$"]
        for j, ax in enumerate(axes.flat):
            ax.legend()
            ax.set_title(titles[j])
            ax.set_xlim(0, 12)
            if j > 1:
                ax.set_xlabel(r"$r$")
        
        # title = get_suptitle(fin, t[i])
        # fig.suptitle(title)

        plt.show()
        plt.close()

        # dq = np.pi * 2 / dx
        # qx = (np.arange(S_rho_AA.shape[1]) - S_rho_AA.shape[1]//2) * dq
        # qy = (np.arange(S_rho_AA.shape[0]) - S_rho_AA.shape[0]//2) * dq
        # qxx, qyy = np.meshgrid(qx, qy)
        # q_module = np.sqrt(qxx ** 2 + qyy ** 2)
        # q_radial = np.arange(S_rho_AA.shape[1]//2) * dq

        n = int(get_Lx(fin) / dx)
        qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        qy = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        q_radial = qx[0:qx.size//2-1]
        qx = np.fft.fftshift(qx)
        qy = np.fft.fftshift(qy)
        qx_ij, qy_ij = np.meshgrid(qx, qy)
        q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)
        S_rho_AA = ori_ave(np.fft.fftshift(S_rho_AA).real, q_radial, q_module)
        S_rho_BB = ori_ave(np.fft.fftshift(S_rho_BB).real, q_radial, q_module)
        S_rho_AB = ori_ave(np.fft.fftshift(S_rho_AB).real, q_radial, q_module)

        S_p_AA = ori_ave(np.fft.fftshift(S_p_AA).real, q_radial, q_module)
        S_p_BB = ori_ave(np.fft.fftshift(S_p_BB).real, q_radial, q_module)
        S_p_AB = ori_ave(np.fft.fftshift(S_p_AB).real, q_radial, q_module)

        S_s_AA = ori_ave(np.fft.fftshift(S_s_AA).real, q_radial, q_module)
        S_s_BB = ori_ave(np.fft.fftshift(S_s_BB).real, q_radial, q_module)
        S_s_AB = ori_ave(np.fft.fftshift(S_s_AB).real, q_radial, q_module)

        S_v_AA = ori_ave(np.fft.fftshift(S_v_AA).real, q_radial, q_module)
        S_v_BB = ori_ave(np.fft.fftshift(S_v_BB).real, q_radial, q_module)
        S_v_AB = ori_ave(np.fft.fftshift(S_v_AB).real, q_radial, q_module)
        fig, axes = plt.subplots(2, 2, figsize=(8, 6.2), constrained_layout=True, sharex=True)
        axes[0, 0].plot(q_radial[1:], S_rho_AA[1:], label="AA")
        axes[0, 0].plot(q_radial[1:], S_rho_BB[1:], label="BB")
        axes[0, 0].plot(q_radial[1:], S_rho_AB[1:], label="AB")

        axes[0, 1].plot(q_radial[1:], S_p_AA[1:], label="AA")
        axes[0, 1].plot(q_radial[1:], S_p_BB[1:], label="BB")
        axes[0, 1].plot(q_radial[1:], S_p_AB[1:], label="AB")

        axes[1, 0].plot(q_radial[1:], S_s_AA[1:], label="AA")
        axes[1, 0].plot(q_radial[1:], S_s_BB[1:], label="BB")
        axes[1, 0].plot(q_radial[1:], S_s_AB[1:], label="AB")

        axes[1, 1].plot(q_radial[1:], S_v_AA[1:], label="AA")
        axes[1, 1].plot(q_radial[1:], S_v_BB[1:], label="BB")
        axes[1, 1].plot(q_radial[1:], S_v_AB[1:], label="AB")
        axes[0, 0].set_xscale("log")

        titles = [r"(a) $\langle \delta \rho_S(\mathbf{q}) \delta \rho_{S'}(-\mathbf{q})\rangle_q$",
                    r"(b) $\langle \mathbf{p}_S(\mathbf{q})\cdot\mathbf{p}_{S'}(-\mathbf{q})\rangle_q$",
                    r"(c) $\langle \delta\bar{v}_S(\mathbf{q})\delta\bar{v}_{S'}(-\mathbf{q})\rangle_q$",
                    r"(d) $\langle \mathbf{v}_S(\mathbf{q})\cdot\mathbf{v}_{S'}(-\mathbf{q})\rangle_q$"]
        for j, ax in enumerate(axes.flat):
            ax.legend()
            ax.set_title(titles[j])
            if j > 1:
                ax.set_xlabel(r"$q$")
        fig.suptitle(title)
        plt.show()
        plt.close()


if __name__ == "__main__":
    dx = 1
    # fin = f"{root_tahmineh}/QS_MPI_noself/L320_r10_10_10_varied_Dt_h/coarse_grain_dx1/L320_320_Dr0.100_Dt0_k0.70_p15_22.5_r10_10_10_e-2.000_J0.500_-0.500_h0.100_1000.npz"
    # beg = 40
    # end = 50
    # fin = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333/coarse_grain_dx1/L640_640_Dr0.100_k0.70_p20_22.5_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    # beg = 20
    # end = 30
    # plot_corr_func(fin, dx, beg, end)

    fin = f"{root_tahmineh}/QS_MPI_noself/L80_Turing_v3_1/coarse_grain_dx1/L80_80_Dr0.100_0.100_k0.70_p40_40_r40_40_40_e-0.200_-2.200_J-2.500_0.400_h0.100_1000.npz"
    beg = 25
    end = 31
    plot_corr_func(fin, dx, beg, end)