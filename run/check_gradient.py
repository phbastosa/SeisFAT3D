import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

path_SPS = "../inputs/geometry/inversion_test_SPS.txt"
path_RPS = "../inputs/geometry/inversion_test_RPS.txt"

nx = 201
ny = 201
nz = 51

dx = 100.0
dy = 100.0
dz = 100.0
 
shot = 1
iteration = 1

# adjoint_grad = pyf.read_binary_volume(nz, nx, ny, f"adjoint_grad_{shot}.bin")
# adjoint_comp = pyf.read_binary_volume(nz, nx, ny, f"adjoint_comp_{shot}.bin")

# gradient = pyf.read_binary_volume(nz-2, nx-2, ny-2, f"gradient_{shot}.bin") 

# adjoint_grad *= 1.0 / np.max(adjoint_grad)
# adjoint_comp *= 1.0 / np.max(adjoint_comp)

# regularization = adjoint_grad / (adjoint_comp + 1e-3)

# scale = np.std(adjoint_grad)

dh = np.array([dx, dy, dz])
slices = np.array([0.6*nz, 0.40*ny, 0.40*nx], dtype = int)

# pyf.plot_model_3D(adjoint_grad, dh, slices, shots = path_SPS, 
#                   nodes = path_RPS, adjx = 0.75, dbar = 1.6,
#                   scale = 2.5, vmin = -scale, vmax = scale)

# pyf.plot_model_3D(adjoint_comp, dh, slices, shots = path_SPS, 
#                   nodes = path_RPS, adjx = 0.75, dbar = 1.6,
#                   scale = 2.5, vmin = -scale, vmax = scale)

# pyf.plot_model_3D(gradient, dh, slices, shots = path_SPS, 
#                   nodes = path_RPS, adjx = 0.75, dbar = 1.6,
#                   scale = 2.5)

# plt.show()