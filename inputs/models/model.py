import numpy as np
import matplotlib.pyplot as plt

nx = 301
ny = 301
nz = 51

dh = 50

true_model = 2500.0*np.ones((nz,nx,ny))

true_model[:40,100:201,100:201] = 2000.0 
true_model[:36,90:211,90:211] = 2000.0
true_model[:32,80:221,80:221] = 2000.0
true_model[:28,70:231,70:231] = 2000.0
true_model[:24,60:241,60:241] = 2000.0
true_model[:20,50:251,50:251] = 2000.0
true_model[:18,:,:] = 2000.0

true_model[:20,100:201,100:201] = 1800.0 
true_model[:18,90:211,90:211] = 1800.0
true_model[:16,80:221,80:221] = 1800.0
true_model[:14,70:231,70:231] = 1800.0
true_model[:12,60:241,60:241] = 1800.0
true_model[:10,:,:] = 1800.0

true_model.flatten("F").astype("float32", order = "F").tofile("true_model_51x301x301_50m.bin")

init_model = np.zeros((nz,nx,ny))

vi = 1800.0
dv = 0.3

z = np.arange(nz) * dh

v = vi + dv*z

for i in range(nz):
    init_model[i] = v[i]

init_model.flatten("F").astype("float32", order = "F").tofile("init_model_51x301x301_50m.bin")

plt.figure(1, figsize = (15,6))
plt.subplot(211)
plt.imshow(true_model[:,:,150], vmin = np.min(true_model), vmax = np.max(true_model))

plt.subplot(212)
plt.imshow(init_model[:,:,150], vmin = np.min(true_model), vmax = np.max(true_model))

plt.tight_layout()
plt.show()