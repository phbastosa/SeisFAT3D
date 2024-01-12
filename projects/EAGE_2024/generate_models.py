import numpy as np

nx = 441
ny = 441
nz = 161

dh = 25

velocity = np.array([2000, 3000, 4000, 5000])
thickness = np.array([1000, 1000, 1000])

accuracy_model = velocity[0] * np.ones((nz, nx, ny)) 

for i in range(len(thickness)):
    accuracy_model[int(np.sum(thickness[:i+1])/dh):] = velocity[i+1]

accuracy_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/accuracyModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
