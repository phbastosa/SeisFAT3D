import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

nrx = 3 
nry = 201

path_SPS = "../inputs/geometry/anisoTest_SPS.txt"
path_RPS = "../inputs/geometry/anisoTest_RPS.txt"
path_XPS = "../inputs/geometry/anisoTest_XPS.txt"

SPS = np.loadtxt(path_SPS, dtype = float, delimiter = ",")
RPS = np.loadtxt(path_RPS, dtype = float, delimiter = ",")
XPS = np.loadtxt(path_XPS, dtype = float, delimiter = ",")

vp = np.array([2000, 3000])
vs = np.array([1180, 1460])
ro = np.array([2250, 2350])
z = np.array([1500])

offset = np.sqrt((SPS[0] - RPS[:,0])**2 + (SPS[1] - RPS[:,1])**2)

directWave = offset / vp[0] 
refracions = pyf.get_analytical_refractions(vp, z, offset)[0]

eikonal_iso = pyf.read_binary_array(nrx*nry,"../outputs/syntheticData/eikonal_iso_nStations603_shot_1.bin")
eikonal_ani = pyf.read_binary_array(nrx*nry,"../outputs/syntheticData/eikonal_ani_nStations603_shot_1.bin")

fig, ax = plt.subplots(ncols = 3, figsize = (15,9))

for n in range(nrx):
    ax[n].plot(directWave[n*nry:n*nry+nry], "o", color = "blue")
    ax[n].plot(refracions[n*nry:n*nry+nry], "o", color = "black")
    
    ax[n].plot(eikonal_iso[n*nry:n*nry+nry], "o", color = "green")
    ax[n].plot(eikonal_ani[n*nry:n*nry+nry], "o", color = "yellow")

    ax[n].set_ylim([0,10])
    ax[n].invert_yaxis()
    
plt.tight_layout()
plt.show()
sys.exit()




