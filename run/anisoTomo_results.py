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
plt.savefig("model_true.png", dpi = 200)
plt.show()

pyf.plot_model_3D(model_init, dh, slices, shots = sps_path, nodes = rps_path, 
                  scale = 2.5, adjx = 0.8, dbar = 1.5, cmap = "jet",
                  cblab = "P wave velocity [km/s]", vmin = 1000, vmax = 4000)
plt.savefig("model_init.png", dpi = 200)
plt.show()

model_pred = pyf.read_binary_volume(nz,nx,ny,"../outputs/models/tomography_iso_final_model_vp_51x201x201.bin")

diff_true = model_true - model_init
diff_pred = model_pred - model_init 

circle_centers = np.array([[8000, 8000],    # SW
                           [8000, 12000],   # NW 
                           [12000, 8000],   # SE
                           [12000, 12000]]) # NE 

pyf.plot_model_3D(diff_true, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_true.png", dpi = 200)
plt.show()

pyf.plot_model_3D(diff_pred, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_pred.png", dpi = 200)
plt.show()

residuo = np.loadtxt("../outputs/residuo/tomography_iso_convergence_5_iterations.txt", dtype = np.float32)

fig, ax = plt.subplots(figsize = (9,4))

ax.plot(100*residuo/np.max(residuo), "--ok")

ax.set_xlabel("Iteration number", fontsize = 15)
ax.set_ylabel("Normalized residuo [%]", fontsize = 15)

ax.set_ylim([0,110])
ax.set_xlim([-0.05,5.05])
fig.tight_layout()
plt.savefig("convergence.png", dpi = 200)
plt.show()
