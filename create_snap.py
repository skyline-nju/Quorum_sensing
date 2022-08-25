import numpy as np
from gsd import hoomd
from read_gsd import read_one_frame


def duplicate(s: hoomd.Snapshot, nx: int, ny: int) -> hoomd.Snapshot:
    N = s.particles.N * nx * ny
    lx = s.configuration.box[0]
    ly = s.configuration.box[1]
    Lx, Ly = lx * nx, ly * ny
    pos = np.zeros((N, 3), dtype=np.float32)
    type_id = np.zeros(N, dtype=np.uint32)
    charge = np.zeros(N, dtype=np.float32)
    for j in range(ny):
        for i in range(nx):
            beg = (j * nx + i) * s.particles.N
            end = beg + s.particles.N
            pos[beg:end, 0] = s.particles.position[:, 0] + lx / 2 + i * lx
            pos[beg:end, 1] = s.particles.position[:, 1] + ly / 2 + j * ly
            pos[beg:end, 2] = s.particles.position[:, 2]
            type_id[beg:end] = s.particles.typeid
            charge[beg:end] = s.particles.charge
    pos[:, 0] -= Lx / 2
    pos[:, 1] -= Ly / 2
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = N
    s2.particles.position = pos
    s2.particles.typeid = type_id
    s2.particles.charge = charge
    s2.particles.types = s.particles.types
    s2.configuration.step = 0
    return s2


def scale(s: hoomd.Snapshot, nx: int, ny: int) -> hoomd.Snapshot:
    N = s.particles.N * nx * ny
    lx = s.configuration.box[0]
    ly = s.configuration.box[1]
    Lx, Ly = lx * nx, ly * ny
    pos = np.zeros((N, 3), dtype=np.float32)
    type_id = np.zeros(N, dtype=np.uint32)
    charge = np.zeros(N, dtype=np.float32)
    for i in range(nx * ny):
        beg = i * s.particles.N
        end = beg + s.particles.N
        pos[beg:end, 0] = s.particles.position[:, 0] * nx
        pos[beg:end, 1] = s.particles.position[:, 1] * ny
        pos[beg:end, 2] = s.particles.position[:, 2]
        type_id[beg:end] = s.particles.typeid
        charge[beg:end] = s.particles.charge
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = N
    s2.particles.position = pos
    s2.particles.typeid = type_id
    s2.particles.charge = charge
    s2.particles.types = s.particles.types
    s2.configuration.step = 0
    return s2


def adjust_density(s: hoomd.Snapshot, phi_A: float, phi_B: float, mode="copy") -> hoomd.Snapshot:
    def add_new_particles(n):
        pos_new = np.zeros((n, 3), dtype=np.float32)
        pos_new[:, 0] = (np.random.rand(n) - 0.5) * Lx
        pos_new[:, 1] = (np.random.rand(n) - 0.5) * Ly
        pos_new[:, 2] = (np.random.rand(n) - 0.5) * np.pi * 2
        return pos_new


    def extend_pos(pos_new, pos_old, mode):
        n_new = pos_new.shape[0]
        n_old = pos_old.shape[0]
        m = n_new // n_old
        n_res = n_new - m * n_old
        end = 0
    
        for i in range(m):
            beg = i * n_old
            end = beg + n_old
            pos_new[beg: end] = pos_old
        
        if n_res > 0:
            if mode == "rand":
                pos_new[end:] = add_new_particles(n_res)
            elif mode == "copy":
                np.random.shuffle(pos_old)
                pos_new[end:] = pos_old[:n_res]
            else:
                print("Error, mode should be copy or rand")
        

    Lx = s.configuration.box[0]
    Ly = s.configuration.box[1]

    n_A = int(Lx * Ly * phi_A)
    n_B = int(Lx * Ly * phi_B)

    pos_A = np.zeros((n_A, 3), dtype=np.float32)
    pos_B = np.zeros((n_B, 3), dtype=np.float32)

    mask_A = s.particles.typeid == 0
    n_A0 = np.sum(mask_A)
    mask_B = s.particles.typeid == 1
    n_B0 = np.sum(mask_B)
    print(s.particles.N, n_A0, n_B0)

    pos_A0 = s.particles.position[mask_A]
    if n_A0 < n_A:
        extend_pos(pos_A, pos_A0, mode)
        print("add %d A particles" % (n_A - n_A0))
    elif n_A0 > n_A:
        np.random.shuffle(pos_A0)
        pos_A = pos_A0[:n_A]
        print("remove %d A particles" % (n_A0 - n_A))
    else:
        pos_A = pos_A0
    
    pos_B0 = s.particles.position[mask_B]
    if n_B0 < n_B:
        extend_pos(pos_B, pos_B0, mode)
        print("add %d B particles" % (n_B - n_B0))
    elif n_B0 > n_B:
        np.random.shuffle(pos_B0)
        pos_B = pos_B0[:n_B]
        print("remove %d B particles" % (n_B0 - n_B))
    else:
        pos_B = pos_B0
    
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = n_A + n_B
    s2.particles.position = np.vstack((pos_A, pos_B))
    s2.particles.typeid = np.zeros(s2.particles.N, dtype=np.float32)
    s2.particles.typeid[n_A:] = 1
    s2.particles.charge = np.zeros(s2.particles.N, dtype=np.float32)
    s2.particles.types = s.particles.types
    s2.configuration.step = 0
    return s2
    

def inverse(x, y, theta, xc, yc):
    x_inv = 2 * xc - x
    y_inv = 2 * yc - y
    theta_inv = theta + np.pi
    theta_inv[theta_inv > np.pi] -= np.pi * 2
    theta_inv[theta_inv < -np.pi] += np.pi * 2
    return x_inv, y_inv, theta_inv


def flip_half(s):
    x, y, theta = s.particles.position.T
    Lx = s.configuration.box[0]
    mask1 = x < 0
    mask2 = x >= 0
    x1, y1, theta1 = x[mask1], y[mask1], theta[mask1]
    x2, y2, theta2 = x[mask2], y[mask2], theta[mask2]
    x1, y1, theta1 = inverse(x1, y1, theta1, -Lx / 4, 0)
    type_id1 = s.particles.typeid[mask1]
    type_id2 = s.particles.typeid[mask2]

    x = np.hstack((x1, x2))
    y = np.hstack((y1, y2))
    theta = np.hstack((theta1, theta2))
    type_id = np.hstack((type_id1, type_id2))
    s.particles.position = np.array([x, y, theta], dtype=np.float32).T
    s.particles.typeid = np.array(type_id, dtype=np.uint32)
    return s


def create_nematic_patten(s, nb=2):
    Lx = s.configuration.box[0]
    s = flip_half(s)
    if nb == 4:
        x = s.particles.position[:, 0]
        mask1 = np.logical_and(x >= -Lx/4, x < 0)
        mask2 = np.logical_and(x >= 0, x <= Lx/4)
        x[mask1] += Lx / 4
        x[mask2] -= Lx / 4
        s.particles.position[:, 0] = x
    s.configuration.step = 0
    return s


if __name__ == "__main__":
    folder = '/scratch03.local/yduan/QS_MPI_noself/L240_80_r40_6_30'
    basename = "L240_80_Dr0.100_k0.70_p3_46_r40_6_30_e-2.000_J2.000_-2.000_20246.gsd"
    fname = f"{folder}/{basename}"
    snap = read_one_frame(fname, -1)
    snap = adjust_density(snap, 4, 46, mode="copy")
    # snap = scale(snap, 2, 2)
    # # snap = create_nematic_patten(snap, nb=4)
    # snap = duplicate(snap, 2, 2)
    snap.configuration.step = 0

    # fout = 'data/snap/L80_80_Dr0.100_k0.70_p38_8_r40_e-2.000_J0.500_-0.500_21000.gsd'
    fout = f"{folder}/L240_80_Dr0.100_k0.70_p4_46_r40_6_30_e-2.000_J2.000_-2.000_20346.gsd"
    f = hoomd.open(name=fout, mode='wb')
    f.append(snap)


