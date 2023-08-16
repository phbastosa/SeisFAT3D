import numpy as np
import matplotlib.pyplot as plt

nx = 501
ny = 501
nz = 51

dh = 50

true_model = 2500.0*np.ones((nz,nx,ny))

true_model[:40,200:301,200:301] = 2000.0 
true_model[:36,190:311,190:311] = 2000.0
true_model[:32,180:321,180:321] = 2000.0
true_model[:28,170:331,170:331] = 2000.0
true_model[:24,160:341,160:341] = 2000.0
true_model[:20,150:351,150:351] = 2000.0
true_model[:18,:,:] = 2000.0

true_model[:20,200:301,200:301] = 1800.0 
true_model[:18,190:311,190:311] = 1800.0
true_model[:16,180:321,180:321] = 1800.0
true_model[:14,170:331,170:331] = 1800.0
true_model[:12,160:341,160:341] = 1800.0
true_model[:10,:,:] = 1800.0

true_model.flatten("F").astype("float32", order = "F").tofile("true_model_51x501x501_50m.bin")

init_model = np.zeros((nz,nx,ny))

vi = 1800.0
dv = 0.3

z = np.arange(nz) * dh

v = vi + dv*z

for i in range(nz):
    init_model[i] = v[i]

init_model.flatten("F").astype("float32", order = "F").tofile("init_model_51x501x501_50m.bin")

plt.figure(1, figsize = (15,6))
plt.subplot(211)
plt.imshow(true_model[:,:,250], vmin = np.min(true_model), vmax = np.max(true_model))

plt.subplot(212)
plt.imshow(init_model[:,:,250], vmin = np.min(true_model), vmax = np.max(true_model))

plt.tight_layout()
plt.show()