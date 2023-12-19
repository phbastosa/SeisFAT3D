import numpy as np
import matplotlib.pyplot as plt

from sys import path

path.append("../src/")

import functions

#---------------------------------------------------------------------

nx = 441
ny = 441
nz = 161

dh = 25.0

shots_file = "../inputs/geometry/xyz_shots_position.txt"
nodes_file = "../inputs/geometry/xyz_nodes_position.txt"
 
shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')

vmin = 2000
vmax = 5000

trueModel = functions.read_binary_volume(nz, nx, ny, f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

subplots = np.array([1, 1], dtype = int)
slices = np.array([int(0.375*nz), int(0.365*nx), int(0.355*ny)], dtype = int) # [xy, zy, zx]
dh = np.array([dh, dh, dh])

functions.plot_model_3D(trueModel, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
plt.savefig(f"true_tomography_model.png", dpi = 200)
plt.clf()

initModel = functions.read_binary_volume(nz, nx, ny, f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dh[0]:.0f}m.bin")

functions.plot_model_3D(initModel, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
plt.savefig(f"init_tomography_model.png", dpi = 200)
plt.clf()

vmin = -500
vmax = 500

functions.plot_model_3D(trueModel-initModel, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
plt.savefig(f"diff_tomography_model.png", dpi = 200)
plt.clf()
