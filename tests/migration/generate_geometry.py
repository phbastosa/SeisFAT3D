import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

parameters = str(sys.argv[1])

nsx = 9
nsy = 3

nrx = 50
nry = 20

ns = nsx*nsy
nr = nrx*nry

sx, sy = np.meshgrid(np.linspace(500, 4500, nsx), 
                     np.linspace(500, 1500, nsy))

rx, ry = np.meshgrid(np.linspace(50, 4950, nrx), 
                     np.linspace(50, 1950, nry))

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
XPS[:, 1] = np.zeros(ns) 
XPS[:, 2] = np.zeros(ns) + nr 

np.savetxt(pyf.catch_parameter(parameters, "SPS"), SPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "RPS"), RPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "XPS"), XPS, fmt = "%.0f", delimiter = ",")
