import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter

nx = 301
ny = 301
nz = 121

dh = 25

true_model = np.zeros((nz,nx,ny))
init_model = np.zeros((nz,nx,ny))

z = np.arange(nz) * dh
vi = 2000.0
dv = 0.675

v = vi + z * dv
v[:40] = 1500

for i in range(nz):
    true_model[i,:,:] = v[i]
    init_model[i,:,:] = v[i]

true_model[60:80,100:200,100:200] = 4000

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/true_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/init_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

shots = np.loadtxt("../inputs/geometry/xyz_shot_positions.txt", delimiter = ",", dtype = np.float32)
nodes = np.loadtxt("../inputs/geometry/xyz_node_positions.txt", delimiter = ",", dtype = np.float32)

gradient = np.fromfile("gradient.bin", dtype = np.float32, count = nx*ny*nz)
gradient = np.reshape(gradient, [nz,nx,ny], order = "F")

new_gradient = np.zeros_like(gradient)

# for i in range(nz):
#     for j in range(nx):
#         for k in range(ny):
#             if i > 0 and i < nz-1 and j > 0 and j < nx-1 and k > 0 and k < ny-1:
#                 dgdx = (gradient[i,j+1,k] - 2*gradient[i,j,k] + gradient[i,j-1,k]) / (dh*dh)
#                 dgdy = (gradient[i,j,k+1] - 2*gradient[i,j,k] + gradient[i,j,k-1]) / (dh*dh)
#                 dgdz = (gradient[i+1,j,k] - 2*gradient[i,j,k] + gradient[i-1,j,k]) / (dh*dh)
                
#                 new_gradient[i,j,k] = dh*dh*np.abs(dgdx + dgdy + dgdz)

# gradient -= new_gradient
gradient = gaussian_filter(gradient, 5.0)

scale = 5.0 * np.std(gradient)

plt.figure(figsize=(15,5))

plt.subplot(131)
plt.imshow(gradient[:,:,150], cmap = "bwr", aspect = "auto", vmin = -scale, vmax = scale)
plt.imshow(true_model[:,:,150], alpha = 0.2, vmin = np.min(true_model), vmax = np.max(true_model))

plt.scatter(shots[:,0]/dh, shots[:,2]/dh)
plt.scatter(nodes[:,0]/dh, nodes[:,2]/dh)

plt.subplot(132)
plt.imshow(gradient[:,150,:], cmap = "bwr", aspect = "auto", vmin = -scale, vmax = scale)
plt.imshow(true_model[:,150,:], alpha = 0.2, vmin = np.min(true_model), vmax = np.max(true_model))

plt.scatter(shots[:,1]/dh, shots[:,2]/dh)
plt.scatter(nodes[:,1]/dh, nodes[:,2]/dh)

plt.subplot(133)
plt.imshow(gradient[70,:,:], cmap = "bwr", aspect = "auto", vmin = -scale, vmax = scale)
plt.imshow(true_model[70,:,:], alpha = 0.2, vmin = np.min(true_model), vmax = np.max(true_model))

plt.scatter(shots[:,0]/dh, shots[:,1]/dh, s = 1)
plt.scatter(nodes[:,0]/dh, nodes[:,1]/dh, s = 1)

plt.tight_layout()
plt.show()
