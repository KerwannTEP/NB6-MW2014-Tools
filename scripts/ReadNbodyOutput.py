import argparse
from scipy.io import FortranFile
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument("-snap_init", help="Initial snapshot available (default: 0)",
                    type=int,default=0)
parser.add_argument("-snap_last", help="Final snapshot available (default: 1000)",
                    type=int,default=1000)

args = parser.parse_args()
snap_init = args.snap_init
snap_last = args.snap_last

# All data is in HU

def read_file(name_old,name_new):
    f = FortranFile(name_old,"r")
    data_int = f.read_record(np.int32)
    data_float = f.read_record(np.float32)

    Npart = data_int[0] # NTOT

    time = data_float[0:1] # TIME [NB]
    pc_per_HU = data_float[2:3] # RBAR [value of one N–body unit in length (pc)]
    Msun_per_HU = data_float[3:4] # ZMBAR [value of one N–body unit in mass (MSun); not necessarily the same value as in the input file]
    r_tidal = data_float[4:5] # RTIDE [tidal radius]
    r_dens = data_float[6:9] # RDENS(1:3) [density center position]
    Myr_per_HU = data_float[10:11] # TSCALE [value of one N–body unit in time (Myr)]
    kms_per_HU = data_float[11:12] # VSTAR [value of one N–body unit in velocity (km/s)]
    r_core = data_float[12:13] #  RC
    n_core = data_float[13:14] # NC [number of stars inside r_core]
    v_core = data_float[14:15] # VC [rms velocity inside r_core]

    # These conversions are in the correct order (they have been tested)

    tab_mass = data_float[0*Npart+20:1*Npart+20].reshape(Npart,1)
    tab_r = data_float[3*Npart+20:6*Npart+20].reshape(Npart,3)
    tab_v = data_float[6*Npart+20:9*Npart+20].reshape(Npart,3)
    tab_pot = data_float[9*Npart+20:10*Npart+20].reshape(Npart,1)

    whead = np.array([time.T[0], pc_per_HU.T[0], Msun_per_HU.T[0], r_tidal.T[0], r_dens.T[0], r_dens.T[1], r_dens.T[2], Myr_per_HU.T[0], kms_per_HU.T[0], r_core.T[0], n_core.T[0], v_core.T[0]])
    w = np.array([tab_mass.T[0], tab_r[:,0].T, tab_r[:,1].T, tab_r[:,2].T, tab_v[:,0].T, tab_v[:,1].T, tab_v[:,2].T, tab_pot.T[0]]).T

    with open(name_new, 'wb') as f_out:
        np.savetxt(f_out, [whead]) # First line of the .txt file
        np.savetxt(f_out, w) # Rest of the .txt file
        

for snap in range(snap_init, snap_last+1, 1):
    name_old = "conf.3_" + str(snap)
    name_new = "output/out_"    + str(snap) + ".txt"
    read_file(name_old,name_new)