import numpy as np

x_max = 2e4
y_max = 2e4
z_max = 5e3

dh = 100.0

nx = int((x_max / dh) + 1)
ny = int((y_max / dh) + 1)
nz = int((z_max / dh) + 1)

init_model = 1500.0 * np.ones((nz, nx, ny))
true_model = 1500.0 * np.ones((nz, nx, ny))

z = np.arange(nz)*dh
dv = 50
vi = 2000

v = vi + z*dv/dh

for i in range(nz):
    true_model[i] = v[i]
    init_model[i] = v[i]

radius = 1000

velocity_variation = np.array([-500, 500, 500, -500])

circle_centers = np.array([[2500, 7500, 7500],
                           [2500, 7500, 11500],
                           [2500, 11500, 7500],
                           [2500, 11500, 11500]])

x, z, y = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh, np.arange(ny)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (y - circle_centers[k,2])**2 + (z - circle_centers[k,0])**2)

    true_model[distance <= radius] += dv

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")