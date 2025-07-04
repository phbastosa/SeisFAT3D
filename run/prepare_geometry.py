import numpy as np

ds = 100.0

offsets = [7000, 8000, 9000]

xc = 10000
yc = 10000

sx = np.array([])
sy = np.array([])

for offset in offsets:

    angles = np.linspace(0.0, 2.0*np.pi, int(2.0*np.pi*offset/ds))

    sx = np.append(sx, offset * np.cos(angles) + xc)
    sy = np.append(sy, offset * np.sin(angles) + yc)


nrx = 19
nry = 19

rx, ry = np.meshgrid(np.linspace(5500, 14500, nrx), np.linspace(5500, 14500, nrx))

# Writing with reciprocity principle

ns = nrx*nry
nr = len(sx)

SPS = np.zeros((ns, 3), dtype = float)
RPS = np.zeros((nr, 3), dtype = float)
XPS = np.zeros((ns, 3), dtype = int)

SPS[:,0] = np.reshape(rx, [ns], order = "F")
SPS[:,1] = np.reshape(ry, [ns], order = "F")
SPS[:,2] = 1000.0

RPS[:,0] = sx
RPS[:,1] = sy
RPS[:,2] = 8.0

XPS[:,0] = np.arange(ns)
XPS[:,1] = np.zeros(ns)
XPS[:,2] = np.zeros(ns) + nr

np.savetxt("../inputs/geometry/anisoTomo_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTomo_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTomo_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")
