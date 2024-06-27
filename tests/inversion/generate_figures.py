import numpy as np
import matplotlib.pyplot as plt

from sys import path
path.append("../src/")
import functions

#-------------------------------------------------------------------------

nx = 201
ny = 201
nz = 51

dh = np.array([100, 100, 100])

true_model = functions.read_binary_volume(nz, nx, ny, f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dh[0]:.0f}m.bin")
init_model = functions.read_binary_volume(nz, nx, ny, f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dh[0]:.0f}m.bin")

shots_file = "../inputs/geometry/xyz_shots_position.txt"
nodes_file = "../inputs/geometry/xyz_nodes_position.txt"
 
slices = np.array([20, 75, 75], dtype = int) # [xy, zy, zx]

functions.plot_model_3D(true_model, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 2000,
                        vmax = 4500,
                        scale = 2.5)

plt.savefig(f"trueModel.png", dpi = 200)
plt.clf()

functions.plot_model_3D(init_model, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 2000,
                        vmax = 4500,
                        scale = 2.5)

plt.savefig(f"initModel.png", dpi = 200)
plt.clf()

functions.plot_model_3D(true_model - init_model, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = -500,
                        vmax = +500,
                        scale = 2.5)

plt.savefig(f"diffModel.png", dpi = 200)
plt.clf()

#-----------------------------------------------

final_model_ls = functions.read_binary_volume(nz, nx, ny, f"../outputs/models/ls_final_model_{nz}x{nx}x{ny}.bin")
final_model_adj = functions.read_binary_volume(nz, nx, ny, f"../outputs/models/adj_final_model_{nz}x{nx}x{ny}.bin")

functions.plot_model_3D(final_model_ls, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 2000,
                        vmax = 4500,
                        scale = 2.5)

plt.savefig(f"final_model_ls.png", dpi = 200)
plt.clf()

functions.plot_model_3D(final_model_adj, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 2000,
                        vmax = 4500,
                        scale = 2.5)

plt.savefig(f"final_model_adj.png", dpi = 200)
plt.clf()

functions.plot_model_3D(final_model_ls - init_model, dh, slices,
                        vmin = -200,
                        vmax = +200,
                        scale = 2.5)

plt.savefig(f"diff_model_ls.png", dpi = 200)
plt.clf()

functions.plot_model_3D(final_model_adj - init_model, dh, slices,
                        vmin = -200,
                        vmax = +200,
                        scale = 2.5)

plt.savefig(f"diff_model_adj.png", dpi = 200)
plt.clf()

# --------------------------------------------------
# model correlation

ls_diff = final_model_ls - true_model
adj_diff = final_model_adj - true_model

rms_error_ls = np.sqrt(np.sum(ls_diff**2)/(nx*ny*nz))
rms_error_adj = np.sqrt(np.sum(adj_diff**2)/(nx*ny*nz))

max_error_ls = np.max(ls_diff)
max_error_adj = np.max(adj_diff)

min_error_ls = np.min(ls_diff)
min_error_adj = np.min(adj_diff)

print(f"|{'-'*74}|")
print(f"| {'Model difference analysis':^30s} | {'RMS ERROR':^11s} | {'MAX ERROR':^11s} | {'MIN ERROR':^11s} |")
print(f"|{'-'*74}|")
print(f"| {'Classical Tomography':^30s} | {f'{rms_error_ls:.4f}':^11s} | {f'{max_error_ls:.5f}':^11s} | {f'{min_error_ls:.4f}':^11s} |")
print(f"| {'Adjoint State Tomography':^30s} | {f'{rms_error_adj:.4f}':^11s} | {f'{max_error_adj:.5f}':^11s} | {f'{min_error_adj:.4f}':^11s} |")
print(f"|{'-'*74}|\n")

shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')  

ns = len(shots)
nr = len(nodes)

ls_data = np.zeros(ns*nr)
adj_data = np.zeros(ns*nr)
obs_data = np.zeros(ns*nr)

for i in range(ns):
    ls_data[i*nr:nr+i*nr] = np.fromfile(f"../outputs/seismograms/ls_fsm_data_nRec400_shot_{i+1}.bin", count = nr, dtype = np.float32)
    adj_data[i*nr:nr+i*nr] = np.fromfile(f"../outputs/seismograms/adj_fsm_data_nRec400_shot_{i+1}.bin", count = nr, dtype = np.float32)
    obs_data[i*nr:nr+i*nr] = np.fromfile(f"../inputs/data/obs_fsm_data_nRec400_shot_{i+1}.bin", count = nr, dtype = np.float32)

rms_error_ls = np.sqrt(np.sum((obs_data - ls_data)**2)/(ns*nr))
rms_error_adj = np.sqrt(np.sum((obs_data - adj_data)**2)/(ns*nr))

max_error_ls = np.max(obs_data - ls_data)
max_error_adj = np.max(obs_data - adj_data)

min_error_ls = np.min(obs_data - ls_data)
min_error_adj = np.min(obs_data - adj_data)

print(f"|{'-'*74}|")
print(f"| {'Data difference analysis':^30s} | {'RMS ERROR':^11s} | {'MAX ERROR':^11s} | {'MIN ERROR':^11s} |")
print(f"|{'-'*74}|")
print(f"| {'Classical Tomography':^30s} | {f'{rms_error_ls:.4f}':^11s} | {f'{max_error_ls:.5f}':^11s} | {f'{min_error_ls:.4f}':^11s} |")
print(f"| {'Adjoint State Tomography':^30s} | {f'{rms_error_adj:.4f}':^11s} | {f'{max_error_adj:.5f}':^11s} | {f'{min_error_adj:.4f}':^11s} |")
print(f"|{'-'*74}|")

# --------------------------------------------------

convergence_ls = np.loadtxt("../outputs/convergence/ls_convergence_5_iterations.txt")
convergence_adj = np.loadtxt("../outputs/convergence/adj_convergence_5_iterations.txt")

plt.figure(21, figsize = (10,4))
plt.plot(convergence_ls, "o--", label = "Least squares approach", color = "orange")
plt.plot(convergence_adj, "o--", label = "Adjoint state approach", color = "green")

plt.title("Convergence curve", fontsize = 18)
plt.xlabel("Iteration number", fontsize = 15)
plt.ylabel(r"$||d^{obs} - d^{cal}||^2_2$", fontsize = 15)

plt.legend(loc = "upper right", fontsize = 12) 
plt.grid(True)
plt.tight_layout()
plt.savefig(f"curve.png", dpi = 200)
