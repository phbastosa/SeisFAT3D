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

model_base = pyf.read_binary_volume(nz,nx,ny,"../outputs/models/base_tomography_iso_final_model_vp_51x201.bin")
model_planeYZ = pyf.read_binary_volume(nz,nx,ny,"../outputs/models/planeYZ_tomography_iso_final_model_vp_51x201.bin")
model_planeXZ = pyf.read_binary_volume(nz,nx,ny,"../outputs/models/planeXZ_tomography_iso_final_model_vp_51x201.bin")
model_planeYZ_rot = pyf.read_binary_volume(nz,nx,ny,"../outputs/models/planeYZ_rot_tomography_iso_final_model_vp_51x201.bin")
model_planeXZ_rot = pyf.read_binary_volume(nz,nx,ny,"../outputs/models/planeXZ_rot_tomography_iso_final_model_vp_51x201.bin")

diff_true = model_true - model_init
diff_base = model_base - model_init 
diff_planeYZ = model_planeYZ - model_init 
diff_planeXZ = model_planeXZ - model_init 
diff_planeYZ_rot = model_planeYZ_rot - model_init 
diff_planeXZ_rot = model_planeXZ_rot - model_init 

circle_centers = np.array([[8000, 8000],    # SW
                           [8000, 12000],   # NW 
                           [12000, 8000],   # SE
                           [12000, 12000]]) # NE 

i = int(0.6*nz)

for center in circle_centers:
    j = int(center[0]/dx)
    k = int(center[1]/dy)

    print(int(diff_base[i,j,k]), int(diff_planeYZ[i,j,k]), int(diff_planeYZ_rot[i,j,k]), int(diff_planeXZ[i,j,k]), int(diff_planeXZ_rot[i,j,k]))


pyf.plot_model_3D(diff_true, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_true.png", dpi = 200)
plt.show()

pyf.plot_model_3D(diff_base, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_base.png", dpi = 200)
plt.show()

pyf.plot_model_3D(diff_planeYZ, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_planeYZ.png", dpi = 200)
plt.show()

pyf.plot_model_3D(diff_planeXZ, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_planeXZ.png", dpi = 200)
plt.show()

pyf.plot_model_3D(diff_planeYZ_rot, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_planeYZ_rot.png", dpi = 200)
plt.show()

pyf.plot_model_3D(diff_planeXZ_rot, dh, slices, scale = 2.6, adjx = 0.82, dbar = 1.5, cmap = "bwr", 
                  cblab = "P wave velocity [km/s]", vmin = -500, vmax = 500)
plt.savefig("diff_planeXZ_rot.png", dpi = 200)
plt.show()

res_base = np.loadtxt("../outputs/residuo/base_tomography_iso_convergence_5_iterations.txt", dtype = np.float32)
res_planeYZ = np.loadtxt("../outputs/residuo/planeYZ_tomography_iso_convergence_5_iterations.txt", dtype = np.float32)
res_planeXZ = np.loadtxt("../outputs/residuo/planeXZ_tomography_iso_convergence_5_iterations.txt", dtype = np.float32)
res_planeYZ_rot = np.loadtxt("../outputs/residuo/planeYZ_rot_tomography_iso_convergence_5_iterations.txt", dtype = np.float32)
res_planeXZ_rot = np.loadtxt("../outputs/residuo/planeXZ_rot_tomography_iso_convergence_5_iterations.txt", dtype = np.float32)

fig, ax = plt.subplots(figsize = (9,4))

ax.plot(100*res_base/np.max(res_base), "--o", label = "Case 0")
ax.plot(100*res_planeYZ/np.max(res_planeYZ), "--o", label = "Case 1")
ax.plot(100*res_planeYZ_rot/np.max(res_planeYZ_rot), "--o", label = "Case 2")
ax.plot(100*res_planeXZ/np.max(res_planeXZ), "--o", label = "Case 3")
ax.plot(100*res_planeXZ_rot/np.max(res_planeXZ_rot), "--o", label = "Case 4")

ax.set_xlabel("Iteration number", fontsize = 15)
ax.set_ylabel("Normalized residuo [%]", fontsize = 15)

ax.set_ylim([0,110])
ax.set_xlim([-0.05,5.05])
ax.legend(loc = "upper right", fontsize = 12)
fig.tight_layout()
plt.savefig("convergence.png", dpi = 200)
plt.show()




