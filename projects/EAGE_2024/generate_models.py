import numpy as np

#---------------------------------------------------------------------------------------------
# Model to show accuracy of travel times
#---------------------------------------------------------------------------------------------

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

#---------------------------------------------------------------------------------------------
# Model to show inaccuracy effects in first arrival tomography
#---------------------------------------------------------------------------------------------

dh = 50

# nx = int((2e4 / dh) + 1)
# ny = int((2e4 / dh) + 1)
# nz = int((5e3 / dh) + 1)

# init_model = np.zeros((nz, nx, ny))
# true_model = np.zeros((nz, nx, ny))

# dv = 50
# vi = 2000

# v = vi + dv*np.arange(nz)

# for i in range(nz):
#     true_model[i] = v[i]
#     init_model[i] = v[i]

# radius = 1000

# velocity_variation = np.array([-500, 500, 500, -500])

# circle_centers = np.array([[2500, 7500, 7500],
#                            [2500, 7500, 11500],
#                            [2500, 11500, 7500],
#                            [2500, 11500, 11500]])

# x, z, y = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh, np.arange(ny)*dh)

# for k, dv in enumerate(velocity_variation):
    
#     distance = np.sqrt((x - circle_centers[k,1])**2 + (y - circle_centers[k,2])**2 + (z - circle_centers[k,0])**2)

#     true_model[distance <= radius] += dv

# true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
# init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")


