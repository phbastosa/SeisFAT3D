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

accuracy_model = 0
#---------------------------------------------------------------------------------------------
# Model to show inaccuracy effects in first arrival tomography
#---------------------------------------------------------------------------------------------

tomography_model = np.zeros((nz, nx, ny))

dv = 15
vi = 2000

v = vi + dv*np.arange(nz)

for i in range(nz):
    tomography_model[i] = v[i]

tomography_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

radius = 1000

velocity_variation = np.array([-500, 500, 500, -500])

circle_centers = np.array([[1500, 4000, 4000],
                           [1500, 4000, 7000],
                           [1500, 7000, 4000],
                           [1500, 7000, 7000]])

x, z, y = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh, np.arange(ny)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (y - circle_centers[k,2])**2 + (z - circle_centers[k,0])**2)

    tomography_model[distance <= radius] += dv

tomography_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")


