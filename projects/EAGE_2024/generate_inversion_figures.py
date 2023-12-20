import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter

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

diffModel = trueModel - initModel

functions.plot_model_3D(diffModel, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
plt.savefig(f"diff_tomography_model.png", dpi = 200)
plt.clf()

fimModel = functions.read_binary_volume(nz, nx, ny, f"../outputs/recovered_models/fim_final_model_{nz}x{nx}x{ny}.bin")

fimDiff = fimModel - initModel

fimDiff = gaussian_filter(fimDiff, 3)

functions.plot_model_3D(fimDiff, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
plt.savefig(f"diff_fim_model.png", dpi = 200)
plt.clf()

fsmModel = functions.read_binary_volume(nz, nx, ny, f"../outputs/recovered_models/fsm_final_model_{nz}x{nx}x{ny}.bin")

fsmDiff = fsmModel - initModel

fsmDiff = gaussian_filter(fsmDiff, 3)

functions.plot_model_3D(fsmDiff, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
plt.savefig(f"diff_fsm_model.png", dpi = 200)
plt.clf()


convergence_fim = np.loadtxt("../outputs/convergence/fim_convergence_10_iterations.txt")
convergence_fsm = np.loadtxt("../outputs/convergence/fsm_convergence_10_iterations.txt")

iterations = np.linspace(0, 10, 11)

plt.figure(1, figsize = (3,6))
plt.plot(convergence_fim, iterations, "-o", color = "orange")
plt.plot(convergence_fsm, iterations, "-o", color = "green")

plt.ylabel("Iteration number", fontsize = 15)
plt.xlabel(r"$||d^{obs} - d^{cal}||^2_2$", fontsize = 15)

plt.yticks(np.arange(11),np.arange(11, dtype = int))

plt.ylim([-0.2,10.2])
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("convergence.png", dpi = 200)




