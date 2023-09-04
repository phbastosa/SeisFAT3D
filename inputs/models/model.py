import numpy as np
import matplotlib.pyplot as plt

nx = 501
ny = 501
nz = 51

dh = 10

vp = 2500.0*np.ones((nz,nx,ny))

vp[:40,200:301,200:301] = 2000.0 
vp[:36,190:311,190:311] = 2000.0
vp[:32,180:321,180:321] = 2000.0
vp[:28,170:331,170:331] = 2000.0
vp[:24,160:341,160:341] = 2000.0
vp[:20,150:351,150:351] = 2000.0
vp[:18,:,:] = 2000.0

vp[:20,200:301,200:301] = 1800.0 
vp[:18,190:311,190:311] = 1800.0
vp[:16,180:321,180:321] = 1800.0
vp[:14,170:331,170:331] = 1800.0
vp[:12,160:341,160:341] = 1800.0
vp[:10,:,:] = 1800.0

vs = vp / 1.7
rho = 310.0*vp**0.25

# vp.flatten("F").astype("float32", order = "F").tofile("vp_model_51x501x501_10m.bin")
# vs.flatten("F").astype("float32", order = "F").tofile("vs_model_51x501x501_10m.bin")
# rho.flatten("F").astype("float32", order = "F").tofile("rho_model_51x501x501_10m.bin")

vmin = np.min(vs)
vmax = np.max(vp)

plt.figure(1, figsize=(15,7))

plt.subplot(311)
plt.imshow(vp[:,:,int(ny/2)], vmin = vmin, vmax = vmax)
plt.colorbar()

plt.subplot(312)
plt.imshow(vs[:,:,int(ny/2)], vmin = vmin, vmax = vmax)
plt.colorbar()

plt.subplot(313)
plt.imshow(rho[:,:,int(ny/2)], vmin = vmin, vmax = vmax)
plt.colorbar()

plt.tight_layout()
plt.show()

init_model = np.zeros_like(vp)

init_model[:25] = vp[0]
init_model[25:] = vp[-1]

delta_vp = vp - init_model 

vp.flatten("F").astype("float32", order = "F").tofile("true_model_51x501x501_10m.bin")
init_model.flatten("F").astype("float32", order = "F").tofile("init_model_51x501x501_10m.bin")

plt.figure(1, figsize = (15,6))
plt.subplot(211)
plt.imshow(vp[:,:,int(ny/2)])

plt.subplot(212)
plt.imshow(init_model[:,:,int(ny/2)])

plt.tight_layout()
plt.show()