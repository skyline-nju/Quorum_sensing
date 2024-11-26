import sys
sys.path.append("../")

from coarse_grain import coarse_grain_all_1D
root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"

if __name__ == "__main__":
    prefix = f"{root_tahmineh}/QS_MPI_noself/surface_instability/L120_40_r10_MA"
    coarse_grain_all_1D(prefix, dx=1)

    prefix = f"{root_tahmineh}/QS_MPI_noself/surface_instability/L120_40_r10_MR"
    coarse_grain_all_1D(prefix, dx=1)