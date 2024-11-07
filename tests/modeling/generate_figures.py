import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

path_SPS = "../inputs/geometry/modeling_test_SPS.txt"
path_RPS = "../inputs/geometry/modeling_test_RPS.txt"

nx = 321
ny = 201
nz = 81 

dx = 25.0
dy = 25.0
dz = 25.0

model_vp = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/modeling_test_vp_model_{nz}x{nx}x{ny}_{dx:.0f}m.bin")
model_vs = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/modeling_test_vs_model_{nz}x{nx}x{ny}_{dx:.0f}m.bin")
model_rho = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/modeling_test_rho_model_{nz}x{nx}x{ny}_{dx:.0f}m.bin")

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model_vp, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.8, dbar = 1.5,
                  cblab = "P wave velocity [km/s]", 
                  vmin = 1000, vmax = 3000)
plt.show()

pyf.plot_model_3D(model_vs, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.8, dbar = 1.5,
                  cblab = "S wave velocity [km/s]",
                  vmin = 1000, vmax = 3000)
plt.show()

pyf.plot_model_3D(model_rho, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.8, dbar = 1.5,
                  cblab = "Density [g/cm³]",
                  vmin = 1000, vmax = 3000)
plt.show()
