import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

path_SPS = "../inputs/geometry/anisoTest_SPS.txt"
path_RPS = "../inputs/geometry/anisoTest_RPS.txt"
path_XPS = "../inputs/geometry/anisoTest_XPS.txt"

SPS = np.loadtxt(path_SPS, dtype = float, delimiter = ",")
RPS = np.loadtxt(path_RPS, dtype = float, delimiter = ",")
XPS = np.loadtxt(path_XPS, dtype = float, delimiter = ",")

nx = 401
ny = 241
nz = 121

dx = 12.5
dy = 12.5
dz = 12.5

model = pyf.read_binary_volume(nz, nx, ny, "../inputs/models/anisoTest_vp.bin")

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.245*ny, 0.247*nx], dtype = int)

pyf.plot_model_3D(model, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, scale = 2.1, adjx = 0.75, 
                  dbar = 1.5, cblab = "P wave velocity [km/s]")
plt.show()

nt = 701
dt = 1e-3
nr = 1500

seismogram = pyf.read_binary_matrix(nt,nr,f"../outputs/syntheticData/elastic_ani_nStations{nr}_nSamples{nt}_shot_1.bin")

scale = 0.1*np.std(seismogram)

fig, ax = plt.subplots(figsize = (12,8))

im = ax.imshow(seismogram[:,int(0.5*nr)-30:int(0.5*nr)+30], aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
cbar = plt.colorbar(im)
cbar.set_label("Amplitude [Pa]", fontsize = 15)

ax.set_xlabel("Receiver index", fontsize = 15)
ax.set_ylabel("Time [s]", fontsize = 15)
ax.set_title("RSG", fontsize = 15)

plt.tight_layout()
plt.show()

