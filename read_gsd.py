from gsd import hoomd, fl
import sys


def get_one_snap(f, i_frame):
    s = hoomd.Snapshot()
    s.configuration.box = f.read_chunk(frame=0, name="configuration/box")

    s.particles.types = f.read_chunk(frame=0, name="particles/types")
    print(s.particles.types)
    try:
        if not isinstance(s.particles.types[0], str):
            s.particles.types = [chr(i) for i in s.particles.types]
    except TypeError:
        s.particles.types = ['A', 'B']
    print(s.particles.types)

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
    s.particles.typeid = f.read_chunk(frame=i_frame,
                                        name="particles/typeid")
    """"
    position = [x, y, theta]
    x \in [-Lx/2, Lx/2)
    y \in [-Ly/2, Ly/2]
    theta \in [-PI, PI]
    """
    s.particles.position = f.read_chunk(frame=i_frame,
                                        name="particles/position")
    return s


def get_nframes(fname):
    with hoomd.open(fname, "rb") as data:
        return len(data)

def read_one_frame(fname, i_frame):
    if i_frame < 0:
        i_frame += get_nframes(fname)
    with fl.open(name=fname, mode="rb") as f:
        return get_one_snap(f, i_frame)


def read_frames(fname, beg=0, end=None, sep=1):
    nframes = get_nframes(fname)
    print("nframes =", nframes)
    if end is None or end > nframes:
        end = nframes
    with fl.open(name=fname, mode="rb") as f:
        for i in range(beg, end, sep):
            snap = get_one_snap(f, i)
            yield snap


if __name__ == "__main__":
    folder = 'QS'
    basename = "L20_40_pA80_pB70_r100_e0.000_a0.800_Dr0.1_Dt0.gsd"
    fname = f"{folder}/{basename}"
    frames = read_frames(fname, beg=0)
    for snap in frames:
        print(snap.configuration.step)



