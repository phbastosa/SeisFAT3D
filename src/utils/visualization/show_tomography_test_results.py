import sys
import numpy as np
import matplotlib.pyplot as plt

from functions import *

filename = sys.argv[1]

x_samples = int(catch_parameter(filename, "x_samples"))
y_samples = int(catch_parameter(filename, "y_samples"))
z_samples = int(catch_parameter(filename, "z_samples"))

x_spacing = float(catch_parameter(filename, "x_spacing"))
y_spacing = float(catch_parameter(filename, "y_spacing"))
z_spacing = float(catch_parameter(filename, "z_spacing"))

true_model_path = catch_parameter(filename, "vp_model_file")

true_model = readBinaryVolume(z_samples, x_samples, y_samples, true_model_path)

shots_file = catch_parameter(filename, "shots_file")
nodes_file = catch_parameter(filename, "nodes_file")

shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([z_samples/2, x_samples/2, y_samples/2], dtype = int) # [xy, zy, zx]
dh = np.array([x_spacing, y_spacing, z_spacing])

# check_geometry(true_model, shots, nodes, dh, slices, subplots)
# plt.savefig("true_model_full.png", dpi = 200)
# plt.show()

# Least squares tomography results

total_iterations = int(catch_parameter(filename, "max_iteration"))
recovered_models_folder = catch_parameter(filename, "estimated_model_folder")
convergence_folder = catch_parameter(filename, "convergence_folder")

ls_model = readBinaryVolume(z_samples, x_samples, y_samples, recovered_models_folder + f"ls_low_model_iteration_2_{z_samples}x{x_samples}x{y_samples}.bin")
adj_model = readBinaryVolume(z_samples, x_samples, y_samples, recovered_models_folder + f"adj_low_model_iteration_2_{z_samples}x{x_samples}x{y_samples}.bin")

check_geometry(true_model - ls_model, shots, nodes, dh, slices, subplots)
plt.savefig(f"ls_diff_model.png", dpi = 200)

check_geometry(true_model - adj_model, shots, nodes, dh, slices, subplots)
plt.savefig(f"adj_diff_model.png", dpi = 200)

ls_diff_max = np.max(true_model - ls_model)
ls_diff_min = np.min(true_model - ls_model)

adj_diff_max = np.max(true_model - adj_model)
adj_diff_min = np.min(true_model - adj_model)

print(f"Max difference between true and ls model = {ls_diff_max}")
print(f"Min difference between true and ls model = {ls_diff_min}")

print(f"Max difference between true and adj model = {adj_diff_max}")
print(f"Min difference between true and adj model = {adj_diff_min}")

# convergence_ls = np.loadtxt(convergence_folder + "ls_low_convergence_2_iterations.txt")
# convergence_adj = np.loadtxt(convergence_folder + "adj_low_convergence_2_iterations.txt")

# n = [0,1,2]

# plt.figure(5, figsize = (10,4))
# plt.plot(n, convergence_ls, "-o", linestyle = "dotted")
# plt.plot(n, convergence_adj, "-o", linestyle = "dotted")
# plt.title("Convergence", fontsize = 18)
# plt.xlabel("Iterations", fontsize = 15)
# plt.ylabel(r"||$\phi(m)$||$^2_2$ = $\sqrt{\sum_n(d_n^{obs} - d_n^{cal})^2}$", fontsize = 15)

# plt.tight_layout()
# plt.savefig("convergence_full.png", dpi = 200)
# plt.show()