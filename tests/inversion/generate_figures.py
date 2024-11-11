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
 
true_model = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/inversion_test_true_model_{nz}x{nx}x{ny}_{dx:.0f}m.bin")
init_model = pyf.read_binary_volume(nz, nx, ny, f"../inputs/models/inversion_test_init_model_{nz}x{nx}x{ny}_{dx:.0f}m.bin")

diff_model = true_model - init_model

dh = np.array([dx, dy, dz])
slices = np.array([0.6*nz, 0.40*ny, 0.40*nx], dtype = int)

pyf.plot_model_3D(true_model, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.75, dbar = 1.6,
                  scale = 2.5, vmin = 1500, vmax = 4000)
plt.savefig("inversion_test_true_model.png", dpi = 300)
plt.show()

pyf.plot_model_3D(init_model, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.75, dbar = 1.6,
                  scale = 2.5, vmin = 1500, vmax = 4000)
plt.savefig("inversion_test_init_model.png", dpi = 300)
plt.show()

pyf.plot_model_3D(diff_model, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.75, dbar = 1.6,
                  scale = 2.5, vmin = -500, vmax = 500, cmap = "bwr")
plt.savefig("inversion_test_diff_model.png", dpi = 300)
plt.show()

least_squares_model = pyf.read_binary_volume(nz, nx, ny, f"../outputs/recoveredModels/inversion_test_least_squares_final_model_{nz}x{nx}x{ny}.bin")
adjoint_state_model = None

least_squares_diff = least_squares_model - init_model
adjoint_state_diff = None

pyf.plot_model_3D(least_squares_model, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.75, dbar = 1.6,
                  scale = 2.5, vmin = 1500, vmax = 4000)
plt.savefig("inversion_test_least_squares_model.png", dpi = 300)
plt.show()

pyf.plot_model_3D(least_squares_diff, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, adjx = 0.75, dbar = 1.6,
                  scale = 2.5, vmin = -500, vmax = 500, cmap = "bwr")
plt.savefig("inversion_test_least_squares_diff.png", dpi = 300)
plt.show()

leastSquaresCurve = np.loadtxt("../outputs/convergence/inversion_test_least_squares_convergence_5_iterations.txt", dtype = float)
adjointStateCurve = None

plt.figure(figsize = (10, 4))

plt.plot(leastSquaresCurve, "--ob", label = "least-squares tomography")

plt.xlabel("Iterations", fontsize = 15)
plt.ylabel(r"$||d^{obs} - d^{cal}||^2_2$", fontsize = 15)
plt.legend(loc = "upper right")

plt.tight_layout()
plt.savefig("inversion_test_convergence.png", dpi = 300)
plt.show()
