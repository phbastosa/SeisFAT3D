import numpy as np
import matplotlib.pyplot as plt

from sys import path, argv
path.append("../src/")
import functions

nx = int(functions.catch_parameter(argv[1], "x_samples"))
ny = int(functions.catch_parameter(argv[1], "y_samples"))
nz = int(functions.catch_parameter(argv[1], "z_samples"))

dx = float(functions.catch_parameter(argv[1], "x_spacing"))
dy = float(functions.catch_parameter(argv[1], "y_spacing"))
dz = float(functions.catch_parameter(argv[1], "z_spacing"))

model = functions.read_binary_volume(nz, nx, ny, functions.catch_parameter(argv[1], "vp_model_file"))

shots_file = functions.catch_parameter(argv[1], "shots_file")
nodes_file = functions.catch_parameter(argv[1], "nodes_file")

slices = np.array([nz/2, ny/2, nx/2], dtype = int)
dh = np.array([dx, dy, dz])

functions.plot_model_3D(model, dh, slices, 
                        shots = shots_file, 
                        nodes = nodes_file,
                        scale = 2.5, 
                        dbar = 1.6)

plt.savefig(f"reference_model.png", dpi = 200)
plt.clf()

model[60:121] = 2300

functions.plot_model_3D(model, dh, slices, 
                        shots = shots_file, 
                        nodes = nodes_file,
                        scale = 2.5, 
                        dbar = 1.6)

plt.savefig(f"initial_hfreq_model.png", dpi = 200)
plt.clf()

model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/init_hfreq_model_{nz}x{nx}x{ny}_{dx:.0f}m.bin")
