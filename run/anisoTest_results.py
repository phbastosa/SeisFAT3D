import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

path_SPS = "../inputs/geometry/anisoTest_SPS.txt"
path_RPS = "../inputs/geometry/anisoTest_RPS.txt"
path_XPS = "../inputs/geometry/anisoTest_XPS.txt"

nx = 201
ny = 101
nz = 51 

dx = 10.0
dy = 10.0
dz = 10.0

model_vp = pyf.read_binary_volume(nz, nx, ny, "../inputs/models/anisoTest_vp.bin")
model_vs = pyf.read_binary_volume(nz, nx, ny, "../inputs/models/anisoTest_vs.bin")
model_ro = pyf.read_binary_volume(nz, nx, ny, "../inputs/models/anisoTest_ro.bin")

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model_vp, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.8, dbar = 1.5,
                  cblab = "P wave velocity [km/s]", 
                  vmin = 1000, vmax = 3000)
plt.savefig("model_homog.png", dpi = 300)
plt.show()

nr = 905

nt = 1001
dt = 1e-3

elastic_iso = pyf.read_binary_matrix(nt,nr,f"../outputs/syntheticData/elastic_iso_nStations{nr}_nSamples{nt}_shot_1.bin")
elastic_vti = pyf.read_binary_matrix(nt,nr,f"../outputs/syntheticData/elastic_vti_nStations{nr}_nSamples{nt}_shot_1.bin")

eikonal_iso = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_iso_nStations{nr}_shot_1.bin")
eikonal_vti = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_vti_nStations{nr}_shot_1.bin")

scale = 5*np.std(elastic_iso)

fig, ax = plt.subplots(num = 10, nrows = 2, figsize = (16,8))

ax[0].imshow(elastic_iso, aspect = "auto", cmap = "Greys", extent = [0, nr, (nt-1)*dt, 0], vmin = -scale, vmax = scale)

ax[0].plot(eikonal_iso, ".g", markersize = 2.0, label = "Eikonal ISO")
ax[0].plot(eikonal_vti, ".r", markersize = 2.0, label = "Eikonal VTI")

ax[0].set_title("ISO", fontsize = 15)
ax[0].set_xlabel("Receiver Index", fontsize = 15)
ax[0].set_ylabel("TWT [s]", fontsize = 15)
ax[0].legend(fontsize = 15)

ax[1].imshow(elastic_vti, aspect = "auto", cmap = "Greys", extent = [0, nr, (nt-1)*dt, 0], vmin = -scale, vmax = scale)

ax[1].plot(eikonal_iso, ".g", markersize = 2.0, label = "Eikonal ISO")
ax[1].plot(eikonal_vti, ".r", markersize = 2.0, label = "Eikonal VTI")

ax[1].set_title("VTI", fontsize = 15)
ax[1].set_xlabel("Receiver Index", fontsize = 15)
ax[1].set_ylabel("TWT [s]", fontsize = 15)
ax[1].legend(fontsize = 15)

fig.tight_layout()
plt.savefig("homog_iso.png", dpi = 300)
plt.show()
