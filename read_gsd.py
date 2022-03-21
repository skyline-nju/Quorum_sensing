# """
# As shown below, errors occur when using hoomd module to read gsd file. Therefore, the gsd file is
# read by fl module instead.

# ---------------------------------------------------------------------------------------------------
# >>> from gsd import hoomd
# >>> fname = "D:/tmp/L40_40_Dr0.100_k0.70_p20_70_r40_e-2.000_J0.020_-0.020_1002.gsd"
# >>> f = hoomd.open(name=fname, mode='rb')
# >>> len(f)
# 297
# >>> s = f[0] 
# Traceback (most recent call last):
#   File "<stdin>", line 1, in <module>
#   File "C:\Users\Yu\miniconda3\envs\py310\lib\site-packages\gsd\hoomd.py", line 1011, in __getitem__
#     return self._read_frame(key)
#   File "C:\Users\Yu\miniconda3\envs\py310\lib\site-packages\gsd\hoomd.py", line 918, in _read_frame
#     tmp = tmp.view(dtype=numpy.dtype((bytes, tmp.shape[1])))
# IndexError: tuple index out of range
# """


from gsd import hoomd, fl
import sys
import numpy as np
import os
import glob


def get_one_snap(f, i_frame):
    s = hoomd.Snapshot()
    s.configuration.box = f.read_chunk(frame=0, name="configuration/box")
    s.particles.types = f.read_chunk(frame=0, name="particles/types")
    try:
        if not isinstance(s.particles.types[0], str):
            s.particles.types = [chr(i) for i in s.particles.types]
    except TypeError:
        s.particles.types = ['A', 'B']

    if f.chunk_exists(frame=i_frame, name="configuration/step"):
        s.configuration.step = f.read_chunk(frame=i_frame,
                                            name="configuration/step")[0]
    else:
        if i_frame == 0:
            s.configuration.step = 0
        else:
            print("Error, cannot find step for frame =", i_frame)
            sys.exit()
    s.particles.N = f.read_chunk(frame=i_frame, name="particles/N")[0]
    s.particles.typeid = f.read_chunk(frame=i_frame, name="particles/typeid")
    """"
    position = [x, y, theta]
    x in [-Lx/2, Lx/2)
    y in [-Ly/2, Ly/2]
    theta in [-PI, PI]
    """
    s.particles.position = f.read_chunk(frame=i_frame,
                                        name="particles/position")

    # 'particles/charge'actually denotes the self-propelled speeds of particles
    if f.chunk_exists(frame=i_frame, name="particles/charge"):
        s.particles.charge = f.read_chunk(frame=i_frame,
                                          name="particles/charge")
    else:
        s.particles.charge = np.zeros(s.particles.N, np.float32)
    return s


def read_one_frame(fname, i_frame):
    with fl.open(name=fname, mode="rb") as f:
        if i_frame < 0:
            i = i_frame + f.nframes
        else:
            i = i_frame
        return get_one_snap(f, i)


def read_frames(fname, beg=0, end=None, sep=1):
    with fl.open(name=fname, mode="rb") as f:
        if end is None or end > f.nframes:
            end = f.nframes
        if beg < 0:
            beg += f.nframes
        for i in range(beg, end, sep):
            snap = get_one_snap(f, i)
            yield snap


def save_last_frames(src_folder):
    dest_folder = src_folder.rstrip("/").split("/")[-1]
    if not os.path.exists(dest_folder):
        os.mkdir(dest_folder)
    
    pat = os.path.join(src_folder, "*.gsd")
    fins = glob.glob(pat)

    for fin in fins:
        snap = read_one_frame(fin, -1)
        fout = os.path.join(dest_folder, os.path.basename(fin))
        with hoomd.open(name=fout, mode="wb") as f:
            f.append(snap)


def get_para(fname):
    def add_para(para, tag, keys, str_list, i):
        status = False
        if tag in str_list[i]:
            s1 = str_list[i].lstrip(tag)
            if s1.replace(".", "").isdigit():
                if isinstance(keys, list):
                    for j, key in enumerate(keys):
                        if j == 0:
                            para[key] = float(s1)
                        else:
                            para[key] = float(str_list[i + j])
                else:
                    para[keys] = float(s1)
                status = True
        return status

    para = {}
    tag_keys = {"L": ["Lx", "Ly"], "Dr": "Dr", "k": "kappa",
                 "p": ["phiA", "phiB"], "r": "rho0", "s": "seed",
                 "e": "eta", "J": ["JAB", "JBA"], "a": "alpha"}
    basename = os.path.basename(fname)
    s = basename.rstrip(".gsd").split("_")

    for i in range(len(s)):
        for tag in tag_keys:
            if add_para(para, tag, tag_keys[tag], s, i):
                break
    if "seed" not in para:
        if s[-1].isdigit():
            para["seed"] = int(s[-1])
    return para
        

if __name__ == "__main__":
    # src_folder = ""
    # save_last_frames(src_folder)

    fname = "data/snap/L80_80_Dr0.100_k0.70_p30_20_r40_e-2.000_J0.300_-0.300_221002.gsd"
    para = get_para(fname)
    print(para)