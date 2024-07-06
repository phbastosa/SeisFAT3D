import numpy as np

from sys import path, argv
path.append("../src/")
import functions

nx = int(functions.catch_parameter(argv[1], "x_samples"))
ny = int(functions.catch_parameter(argv[1], "y_samples"))
nz = int(functions.catch_parameter(argv[1], "z_samples"))

dz = float(functions.catch_parameter(argv[1], "z_spacing"))

model_vp = np.zeros((nz,nx,ny)) + 1500
model_vs = np.zeros((nz,nx,ny))
model_rho = np.zeros((nz,ny,nx)) + 1000
        
model_vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_vp_{nz}x{nx}x{ny}_{dz:.0f}m.bin")
model_vs.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_vs_{nz}x{nx}x{ny}_{dz:.0f}m.bin")
model_rho.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_rho_{nz}x{nx}x{ny}_{dz:.0f}m.bin")
