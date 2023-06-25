import numpy as np
import matplotlib.pyplot as plt

nx = 301
ny = 301
nz = 121

dh = 25

true_model = np.zeros((nz,nx,ny))
init_model = np.zeros((nz,nx,ny))

z = np.arange(nz) * dh
vi = 2000.0
dv = 0.675

v = vi + z * dv

for i in range(nz):
    true_model[i,:,:] = v[i]
    init_model[i,:,:] = v[i]

true_model[40:80,100:200,100:200] = 2500

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/true_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/init_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
