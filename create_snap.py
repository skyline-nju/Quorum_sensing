import numpy as np
from gsd import hoomd, fl
import matplotlib.pyplot as plt


def read_snap(fname, i_frame):
    with fl.open(name=fname, mode="rb") as f:
        s = hoomd.Snapshot()
        s.configuration.box = f.read_chunk(frame=0, name="configuration/box")
        s.particles.types = f.read_chunk(frame=0, name="particles/types")
        s.configuration.step = f.read_chunk(frame=i_frame,
                                            name="configuration/step")[0]
        s.particles.N = f.read_chunk(frame=i_frame, name="particles/N")[0]
        s.particles.typeid = f.read_chunk(frame=i_frame,
                                          name="particles/typeid")
        # position = [x, y, theta]
        # x \in [-Lx/2, Lx/2)
        # y \in [-Ly/2, Ly/2]
        # theta \in [-PI, PI]
        s.particles.position = f.read_chunk(frame=i_frame,
                                            name="particles/position")
    return s


def duplicate(s, nx, ny):
    N = s.particles.N * nx * ny
    lx = s.configuration.box[0]
    ly = s.configuration.box[1]
    Lx, Ly = lx * nx, ly * ny
    pos = np.zeros((N, 3), dtype=np.float32)
    type_id = np.zeros(N, dtype=np.uint32)
    for j in range(ny):
        for i in range(nx):
            beg = (j * nx + i) * s.particles.N
            end = beg + s.particles.N
            pos[beg:end, 0] = s.particles.position[:, 0] + lx / 2 + i * lx
            pos[beg:end, 1] = s.particles.position[:, 1] + ly / 2 + j * ly
            pos[beg:end, 2] = s.particles.position[:, 2]
            type_id[beg:end] = s.particles.typeid
    pos[:, 0] -= Lx / 2
    pos[:, 1] -= Ly / 2
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.types = ['A', 'B']
    s2.particles.N = N
    s2.particles.position = pos
    s2.particles.typeid = type_id
    s2.configuration.step = 0
    return s2


def inverse(x, y, theta, xc, yc):
    x_inv = 2 * xc - x
    y_inv = 2 * yc - y
    theta_inv = theta + np.pi
    theta_inv[theta_inv > np.pi] -= np.pi * 2
    theta_inv[theta_inv < -np.pi] += np.pi * 2
    return x_inv, y_inv, theta_inv


def create_polar_patten(s):
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
    s.particles.types = ['A', 'B']
    return s


if __name__ == "__main__":
    fname = "D:/data/QS/20_5_200_100_0.00_2.00_1.0_0.02_0.gsd"
    snap = read_snap(fname, 1832)

    # snap2 = duplicate(snap, 1, 4)

    snap = create_polar_patten(snap)

    f = hoomd.open(name='b.gsd', mode='wb')
    f.append(snap)
