import numpy as np

dh = np.array([100, 50, 25])

nx = np.array(22000.0 / dh + 1, dtype = int)
ny = np.array(22000.0 / dh + 1, dtype = int)
nz = np.array( 5000.0 / dh + 1, dtype = int)

interfaces = np.array([1000, 2500, 4500])
velocities = np.array([2000, 3000, 4000])

for i, spacing in enumerate(dh):
    model = 1500.0*np.ones((nz[i], nx[i], ny[i]))
    
    for j, depth in enumerate(interfaces):
        model[int(depth/spacing):] = velocities[j]

    model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/testModel_{int(nz[i])}x{int(nx[i])}x{int(ny[i])}_{spacing:.0f}m.bin")
