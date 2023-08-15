import numpy as np
import matplotlib.pyplot as plt

nx = 301
ny = 301
nz = 51

model = 2500.0*np.ones((nz,nx,ny))

model[:40,100:201,100:201] = 2000.0 
model[:36,90:211,90:211] = 2000.0
model[:32,80:221,80:221] = 2000.0
model[:28,70:231,70:231] = 2000.0
model[:24,60:241,60:241] = 2000.0
model[:20,50:251,50:251] = 2000.0
model[:18,:,:] = 2000.0

model[:20,100:201,100:201] = 1800.0 
model[:18,90:211,90:211] = 1800.0
model[:16,80:221,80:221] = 1800.0
model[:14,70:231,70:231] = 1800.0
model[:12,60:241,60:241] = 1800.0
model[:10,:,:] = 1800.0

model.flatten("F").astype("float32", order = "F").tofile("true_model_51x301x301_50m.bin")

model = 2500.0*np.ones((nz,nx,ny))
model[:25] = 2000.0

model.flatten("F").astype("float32", order = "F").tofile("init_model_51x301x301_50m.bin")

plt.imshow(model[:,:,150])
plt.show()