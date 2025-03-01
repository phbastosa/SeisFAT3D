import numpy as np
import matplotlib.pyplot as plt

nx = 201
ny = 201
nz = 201

P = np.fromfile("P.bin", dtype = np.float32, count = nx*ny*nz).reshape([nz,nx,ny], order = "F")

P = P[50:-50,50:-50,50:-50]

nx = 105
ny = 105
nz = 105

T = np.fromfile("T.bin", dtype = np.float32, count = nx*ny*nz).reshape([nz,nx,ny], order = "F")

T = T[2:-2,2:-2,2:-2]

nx = 101
ny = 101
nz = 101

plt.imshow(P[:,:,int(0.5*ny)])
plt.contour(T[:,:,int(0.5*ny)], levels = [0.2])
plt.show()
