import numpy as np

nsx = 2
nsy = 2

nrx = 157 
nry = 4

ns = nsx*nsy
nr = nrx*nry

sx, sy = np.meshgrid(np.linspace(100, 7900, nsx), np.linspace(100, 4900, nsy))

rx, ry = np.meshgrid(np.linspace(100, 7900, nrx), np.linspace(1000, 4000, nry))

SPS = np.zeros((ns, 3), dtype = float)
RPS = np.zeros((nr, 3), dtype = float)
XPS = np.zeros((ns, 3), dtype = int)

SPS[:,0] = np.reshape(sx, [ns], order = "F")
SPS[:,1] = np.reshape(sy, [ns], order = "F")
SPS[:,2] = np.zeros(ns)

RPS[:,0] = np.reshape(rx, [nr], order = "C")
RPS[:,1] = np.reshape(ry, [nr], order = "C")
RPS[:,2] = np.zeros(nr)

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.arange(ns)*nrx 
XPS[:, 2] = np.arange(ns)*nrx + nrx 

np.savetxt("../inputs/geometry/modeling_test_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/modeling_test_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/modeling_test_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")
