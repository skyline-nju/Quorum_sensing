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
    folder = '/scratch03.local/yduan/QS5/polar_bands'
    basename = "L160_5_Dr0.100_k0.70_p40_40_r40_e-0.900_J0.500_-0.500_411200.gsd"
    fname = f"{folder}/{basename}"
    snap = read_one_frame(fname, -1)
    # snap = create_nematic_patten(snap, nb=4)
    snap = duplicate(snap, 2, 1)
    # snap.configuration.step = 0

    fout = 'data/snap/L320_5_Dr0.100_k0.70_p40_40_r40_e-0.900_J0.500_-0.500_811200.gsd'
    # fout = f"{folder}/L40_5_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.500_-0.500_41200.gsd"
    f = hoomd.open(name=fout, mode='wb')
    f.append(snap)
