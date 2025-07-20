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
    NTOT, MODEL, NRUN, NK = data_int
    
    # NTOT = N + NPAIRS (Number of single particles plus centres of masses of regularized (KS) pairs)
    
    # read mixed block as raw bytes:
    # To avoid issues with the NAME(J) entry
    # This is because we switched from float to doubles in output.F for everything but NAME(J)
    raw = f.read_record(np.uint8)  # raw bytes
    nb = 8
    
    # Single entries
    AD = np.frombuffer(raw[0:nb*(NK)], dtype=np.float64)
    
    time = AD[0:1] # TIME [NB]
    pc_per_HU = AD[2:3] # RBAR [value of one N–body unit in length (pc)]
    Msun_per_HU = AD[3:4] # ZMBAR [value of one N–body unit in mass (MSun); not necessarily the same value as in the input file]
    r_tidal = AD[4:5] # RTIDE [tidal radius]
    r_dens = AD[6:9] # RDENS(1:3) [density center position]
    Myr_per_HU = AD[10:11] # TSCALE [value of one N–body unit in time (Myr)]
    kms_per_HU = AD[11:12] # VSTAR [value of one N–body unit in velocity (km/s)]
    r_core = AD[12:13] #  RC
    n_core = AD[13:14] # NC [number of stars inside r_core]
    v_core = AD[14:15] # VC [rms velocity inside r_core]
    r_cluster = AD[20:23] # RG(1:3)
    v_cluster = AD[23:26] # VG(1:3)
    
    # Following arrays arrays are of size N + NPAIRS
    # The last line is the c.m. : a ghost particle
    # We only save the first N lines
    
    NPAIRS = int(AD[1])
    N = NTOT - NPAIRS

    tab_mass = np.frombuffer(raw[nb*(NK):nb*(NK+NTOT)], dtype=np.float32).reshape(NTOT,1)[0:N,:]
    tab_r = np.frombuffer(raw[nb*(NK+3*NTOT):nb*(NK+6*NTOT)], dtype=np.float32).reshape(NTOT,3)[0:N,:]
    tab_v = np.frombuffer(raw[nb*(NK+6*NTOT):nb*(NK+9*NTOT)], dtype=np.float32).reshape(NTOT,3)[0:N,:]
    tab_pot = np.frombuffer(raw[nb*(NK+9*NTOT):nb*(NK+10*NTOT)], dtype=np.float32).reshape(NTOT,1)[0:N,:]

    whead = np.array([time.T[0], pc_per_HU.T[0], Msun_per_HU.T[0], r_tidal.T[0], r_dens.T[0], r_dens.T[1], r_dens.T[2], Myr_per_HU.T[0], kms_per_HU.T[0], r_core.T[0], n_core.T[0], v_core.T[0], r_cluster.T[0], r_cluster.T[1], r_cluster.T[2], v_cluster.T[0], v_cluster.T[1], v_cluster.T[2]])
    w = np.array([tab_mass.T[0], tab_r[:,0].T, tab_r[:,1].T, tab_r[:,2].T, tab_v[:,0].T, tab_v[:,1].T, tab_v[:,2].T, tab_pot.T[0]]).T

    with open(name_new, 'wb') as f_out:
        np.savetxt(f_out, [whead]) # First line of the .txt file
        np.savetxt(f_out, w) # Rest of the .txt file
        

for snap in range(snap_init, snap_last+1, 1):
    name_old = "conf.3_" + str(snap)
    name_new = "output/out_"    + str(snap) + ".txt"
    read_file(name_old,name_new)