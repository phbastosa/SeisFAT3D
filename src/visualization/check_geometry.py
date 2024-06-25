import numpy as np
import matplotlib.pyplot as plt

from sys import path,argv
path.append("../src/")
import functions 

nx = int(functions.catch_parameter(argv[1], "x_samples"))
ny = int(functions.catch_parameter(argv[1], "y_samples"))
nz = int(functions.catch_parameter(argv[1], "z_samples"))

dx = float(functions.catch_parameter(argv[1], "x_spacing"))
dy = float(functions.catch_parameter(argv[1], "y_spacing"))
dz = float(functions.catch_parameter(argv[1], "z_spacing"))

model_file = functions.catch_parameter(argv[1], "vp_model_file")

model = functions.read_binary_volume(nz, nx, ny, model_file)

shots_file = functions.catch_parameter(argv[1], "shots_file")
nodes_file = functions.catch_parameter(argv[1], "nodes_file")

slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

dh = np.array([dx, dy, dz])

functions.plot_model_3D(model, dh, slices, 
                        shots = shots_file, 
                        nodes = nodes_file,
                        scale = 2.5, 
                        dbar = 1.6)

plt.savefig(f"current_vp_model.png", dpi = 200)
print("\nFile current_vp_model.png was succesfully written.")
