import struct
import numpy as np
import matplotlib.pyplot as plt


def read_one_file(fin, Lx, Ly, cg_l, cg_dt):
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()
        nx, ny = Lx // cg_l, Ly // cg_l
        n = nx * ny
        framesize = n * 4
        nframe = filesize // framesize
        print("find", nframe, "frames")
        if nframe * framesize != filesize:
            print("Warning, filesize for", fin, "is", filesize,
                  ", but should be to be", nframe * framesize)
        f.seek(0)
        file_end = filesize
        frame_idx = 0
        while f.tell() < file_end:
            buf = f.read(framesize)
            data = np.array(struct.unpack("%dH" % (n * 2), buf))
            rho_A, rho_B = data.reshape(2, ny, nx) / cg_l**2
            frame_idx += 1
            t = frame_idx * cg_dt
            yield t, rho_A, rho_B


if __name__ == "__main__":
    root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_r10_10_10_h0.333/coarse_grain_dx4"


    fin = f"{folder}/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2040_t1000_0.bin"

    frames = read_one_file(fin, 1280, 1280, 4, 1000)
    for (t, rho_A, rho_B) in frames:
        plt.imshow(rho_B, origin="lower")
        plt.show()
        plt.close()