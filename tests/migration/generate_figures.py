import sys; sys.path.append("../src/")

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import functions as pyf

path_SPS = "../inputs/geometry/migration_test_SPS.txt"
path_RPS = "../inputs/geometry/migration_test_RPS.txt"

nx = 241
ny = 161
nz = 41 

dx = 12.5
dy = 12.5
dz = 12.5

model_vp = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/migration_test_vp_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
model_vs = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/migration_test_vs_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
model_rho = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/migration_test_rho_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.37*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model_vp, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.85, dbar = 1.7,
                  scale = 4.3, cblab = "P wave velocity [km/s]", 
                  vmin = 1000, vmax = 2000)
plt.savefig("migration_test_vp.png", dpi = 300)

pyf.plot_model_3D(model_vs, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.85, dbar = 1.7,
                  scale = 4.3,cblab = "S wave velocity [km/s]",
                  vmin = 1000, vmax = 2000)
plt.savefig("migration_test_vs.png", dpi = 300)

pyf.plot_model_3D(model_rho, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.85, dbar = 1.7,
                  scale = 4.3, cblab = "Density [g/cm³]",
                  vmin = 1000, vmax = 2000)
plt.savefig("migration_test_rho.png", dpi = 300)

image = pyf.read_binary_volume(nz, nx, ny, f"../outputs/migratedImages/kirchhoff_result_{nz}x{nx}x{ny}.bin")

scale = 5.0*np.std(image)

slices = np.array([0.5*nz, 0.37*ny, 0.5*nx], dtype = int)

image = sc.ndimage.laplace(image)

pyf.plot_model_3D(image, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.85, dbar = 1.7,
                  scale = 4.3, cblab = "Amplitude",
                  vmin = -scale*1.5, vmax = scale*1.5, cmap = "seismic")
plt.savefig("migration_test_result.png", dpi = 300)
