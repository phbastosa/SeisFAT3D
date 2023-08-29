import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter

nz = 51 
nx = 501
ny = 501

data = np.fromfile("../outputs/gradient/gradient_iteration_0_51x501x501.bin", dtype = np.float32, count = nx*ny*nz)
gradient = np.reshape(data, [nz,nx,ny], order = "F")

# gradient = gaussian_filter(gradient, 5.0)

sx = 50
sy = 50

# gradient[:,:,sy] = 0.5*(gradient[:,:,sy+1] + gradient[:,:,sy-1])
# gradient[:,sx,:] = 0.5*(gradient[:,sx+1,:] + gradient[:,sx-1,:])

plt.figure(1, figsize=(5,5))

scale = np.std(gradient[25,:,:])

plt.imshow(gradient[25,:,:], aspect = "auto")

plt.tight_layout()
plt.show()
