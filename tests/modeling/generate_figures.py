import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

path_SPS = "../inputs/geometry/modeling_test_SPS.txt"
path_RPS = "../inputs/geometry/modeling_test_RPS.txt"
path_XPS = "../inputs/geometry/modeling_test_XPS.txt"

nx = 201
ny = 201
nz = 201 

dx = 50.0
dy = 50.0
dz = 50.0

model_vp = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/modeling_test_vp.bin")

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model_vp, dh, slices, shots = path_SPS, scale = 0.4, 
                  nodes = path_RPS, adjx = 0.5, dbar = 1.25,
                  cblab = "P wave velocity [km/s]", 
                  vmin = 1600, vmax = 2000)
plt.savefig("modeling_test_vp.png", dpi = 300)

SPS = np.loadtxt(path_SPS, dtype = float, delimiter = ",")
RPS = np.loadtxt(path_RPS, dtype = float, delimiter = ",")
XPS = np.loadtxt(path_XPS, dtype = int, delimiter = ",")

nr = len(RPS)

eikonal_iso = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_iso_nStations{nr}_shot_1.bin")
eikonal_ani = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_ani_nStations{nr}_shot_1.bin")

offset = np.sqrt((SPS[0] - RPS[:,0])**2 + (SPS[1] - RPS[:,1])**2 + (SPS[2] - RPS[:,2])**2)

analyticalT = offset / model_vp[0,0,0] 

fig, ax = plt.subplots(figsize = (16,6))
  
ax.plot(eikonal_iso, ".", label = "Eikonal Isotropic")
ax.plot(eikonal_ani, ".", label = r"Eikonal Anisotropic $\epsilon = 0.1$")
ax.plot(analyticalT, ".", label = "Analytical isotropic")

ax.set_ylabel("Time [s]", fontsize = 15)
ax.set_xlabel("Trace index", fontsize = 15)
ax.legend(loc = "upper right", fontsize = 15)

ax.invert_yaxis()
fig.tight_layout()
plt.savefig("modeling_test_data.png", dpi = 300)
