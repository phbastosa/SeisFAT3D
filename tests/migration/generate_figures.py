import sys; sys.path.append("../src/")

import numpy as np
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
slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model_vp, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.85, dbar = 1.7,
                  scale = 4.3, cblab = "P wave velocity [km/s]", 
                  vmin = 1500, vmax = 2500)
plt.savefig("modeling_test_vp.png", dpi = 300)

pyf.plot_model_3D(model_vs, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.85, dbar = 1.7,
                  scale = 4.3,cblab = "S wave velocity [km/s]",
                  vmin = 1500, vmax = 2500)
plt.savefig("modeling_test_vs.png", dpi = 300)

pyf.plot_model_3D(model_rho, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.85, dbar = 1.7,
                  scale = 4.3, cblab = "Density [g/cm³]",
                  vmin = 1500, vmax = 2500)
plt.savefig("modeling_test_rho.png", dpi = 300)
plt.show()
# nt = 5001
# dt = 1e-3

# ns = 4
# nr = 157

# fig, ax = plt.subplots(ncols = 4, figsize = (16,6))

# for i in range(ns):
    
#     eikonal = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_iso_nStations157_shot_{i+1}.bin")
#     elastic = pyf.read_binary_matrix(nt, nr, f"../outputs/syntheticData/elastic_iso_nStations157_nSamples5001_shot_{i+1}.bin")

#     scale = 0.9*np.std(elastic)

#     ax[i].imshow(elastic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
#     ax[i].plot(eikonal / dt, "--")

# plt.show()
