import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

parameters = str(sys.argv[1])

sps_path = pyf.catch_parameter(parameters, "SPS") 
rps_path = pyf.catch_parameter(parameters, "RPS") 
xps_path = pyf.catch_parameter(parameters, "XPS") 

nx = int(pyf.catch_parameter(parameters, "x_samples"))
ny = int(pyf.catch_parameter(parameters, "y_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples"))

dx = float(pyf.catch_parameter(parameters, "x_spacing"))
dy = float(pyf.catch_parameter(parameters, "y_spacing"))
dz = float(pyf.catch_parameter(parameters, "z_spacing"))

model_vp = pyf.read_binary_volume(nz, nx, ny, pyf.catch_parameter(parameters, "vp_model_file"))

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model_vp, dh, slices, shots = sps_path, nodes = rps_path, scale = 14, 
                 adjx = 0.8, dbar = 1.5, cmap = "Greys", cblab = "Vp [m/s]")
plt.savefig("migration_test_setup.png", dpi = 200)
plt.show()

