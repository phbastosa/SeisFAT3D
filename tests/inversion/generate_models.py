import numpy as np

x_max = 10000
y_max = 10000
z_max = 6000

dh = 50.0

nx = int((x_max / dh) + 1)
ny = int((y_max / dh) + 1)
nz = int((z_max / dh) + 1)

init_model = 1500.0 * np.ones((nz, nx, ny))
true_model = 1500.0 * np.ones((nz, nx, ny))

radius = 1000

velocity_variation = np.array([-500, 500, 500, -500])

circle_centers = np.array([[3000, 3500, 3500],
                           [3000, 3500, 6500],
                           [3000, 6500, 3500],
                           [3000, 6500, 6500]])

x, z, y = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh, np.arange(ny)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (y - circle_centers[k,2])**2 + (z - circle_centers[k,0])**2)

    true_model[distance <= radius] += dv

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")