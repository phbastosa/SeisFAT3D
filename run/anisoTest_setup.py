import numpy as np

nx = 201
ny = 101
nz = 51

dh = 10.0

Vp = np.zeros((nz,nx,ny)) + 2000.0
Vs = np.zeros((nz,nx,ny)) + 1100.0
Ro = np.zeros((nz,nx,ny)) + 2600.0

E = np.zeros((nz,nx,ny))
D = np.zeros((nz,nx,ny))
G = np.zeros((nz,nx,ny))

Vp[25:] = 2500.0
Vs[25:] = 1450.0
Ro[25:] = 2630.0

E[25:] = 0.0
D[25:] = 0.0
G[25:] = -0.2

Vp.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_vp.bin")
Vs.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_vs.bin")
Ro.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_ro.bin")

E.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_E.bin")
D.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_D.bin")
G.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_G.bin")

###############

nsx = 1
nsy = 1

nrx = 181 
nry = 5

ns = nsx*nsy
nr = nrx*nry

sx, sy = 50, 500 
rx, ry = np.meshgrid(np.linspace(100, 1900, nrx), np.linspace(100, 900, nry))

SPS = np.zeros((ns, 3), dtype = float)
RPS = np.zeros((nr, 3), dtype = float)
XPS = np.zeros((ns, 3), dtype = int)

SPS[:,0] = np.reshape(sx, [ns], order = "F")
SPS[:,1] = np.reshape(sy, [ns], order = "F")
SPS[:,2] = np.zeros(ns) + dh

RPS[:,0] = np.reshape(rx, [nr], order = "C")
RPS[:,1] = np.reshape(ry, [nr], order = "C")
RPS[:,2] = np.zeros(nr) + dh*(nz-1) - dh

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.zeros(ns) 
XPS[:, 2] = np.zeros(ns) + nr 

np.savetxt("../inputs/geometry/anisoTest_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTest_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTest_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")
