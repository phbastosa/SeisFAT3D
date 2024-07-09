import numpy as np

nx = 201
ny = 201
nz = 201

dz = 10.0

model_vp = np.zeros((nz,nx,ny)) + 1500
model_vs = np.zeros((nz,nx,ny))
model_rho = np.zeros((nz,ny,nx)) + 1000
        
model_vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_vp_{nz}x{nx}x{ny}_{dz:.0f}m.bin")
model_vs.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_vs_{nz}x{nx}x{ny}_{dz:.0f}m.bin")
model_rho.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_rho_{nz}x{nx}x{ny}_{dz:.0f}m.bin")
