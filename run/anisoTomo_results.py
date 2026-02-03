import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf
import matplotlib.pyplot as plt

parameters = str(sys.argv[1])

sps_path = pyf.catch_parameter(parameters,"SPS") 
rps_path = pyf.catch_parameter(parameters,"RPS") 

nx = int(pyf.catch_parameter(parameters, "x_samples"))
ny = int(pyf.catch_parameter(parameters, "y_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples")) 

dx = float(pyf.catch_parameter(parameters, "x_spacing"))
dy = float(pyf.catch_parameter(parameters, "y_spacing"))
dz = float(pyf.catch_parameter(parameters, "z_spacing"))

model_true = np.fromfile(f"../inputs/models/anisoTomo_vp_true.bin", count = nx*ny*nz, dtype = np.float32).reshape([nz,nx,ny], order = "F")
model_init = np.fromfile(f"../inputs/models/anisoTomo_vp_init.bin", count = nx*ny*nz, dtype = np.float32).reshape([nz,nx,ny], order = "F")

SPS = np.loadtxt("../inputs/geometry/anisoTomo_SPS.txt", delimiter = ",", dtype = float)
RPS = np.loadtxt("../inputs/geometry/anisoTomo_RPS.txt", delimiter = ",", dtype = float)

dh = np.array([dz, dy, dx])
slices = np.array([0.6*nz, 0.4*ny, 0.4*nx], dtype = int)

pyf.plot_model_3D(model_true, dh, slices, shots = sps_path, nodes = rps_path, 
                  scale = 2.5, adjx = 0.8, dbar = 1.5, cmap = "jet",
                  cblab = "P wave velocity [km/s]", vmin = 1000, vmax = 4000)
plt.show()

pyf.plot_model_3D(model_init, dh, slices, shots = sps_path, nodes = rps_path, 
                  scale = 2.5, adjx = 0.8, dbar = 1.5, cmap = "jet",
                  cblab = "P wave velocity [km/s]", vmin = 1000, vmax = 4000)
plt.show()

sId = 200
nr = 1506

dobs_iso = pyf.read_binary_array(nr, f"../inputs/data/obs_eikonal_iso_nStations1506_shot_{sId}.bin")
dobs_ani = pyf.read_binary_array(nr, f"../inputs/data/obs_eikonal_ani_nStations1506_shot_{sId}.bin")

fig, ax = plt.subplots()

ax.plot(dobs_iso, "ok")
ax.plot(dobs_ani, "or")

fig.tight_layout()
plt.show()
