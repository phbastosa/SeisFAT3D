import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

nt = 1001
dt = 1e-3

ns = 164
nr = 200

RPS = np.loadtxt("../inputs/geometry/migration_test_RPS.txt", delimiter = ",", dtype = float)
SPS = np.loadtxt("../inputs/geometry/migration_test_SPS.txt", delimiter = ",", dtype = float)
XPS = np.loadtxt("../inputs/geometry/migration_test_XPS.txt", delimiter = ",", dtype = int)

for s in range(ns):

    sx = SPS[s,0]
    sy = SPS[s,1]
    sz = SPS[s,2]

    rx = RPS[XPS[s,1]:XPS[s,2], 0]
    ry = RPS[XPS[s,1]:XPS[s,2], 1]
    rz = RPS[XPS[s,1]:XPS[s,2], 2]

    distance = np.sqrt((sx - rx)**2 + (sy - ry)**2 + (sz - rz)**2)

    time = distance / 1200.0

    ts = np.array(time/dt + 0.1/dt, dtype = int)

    elastic = pyf.read_binary_matrix(nt, nr, f"../inputs/data/migration_test_elastic_iso_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")
    
    for r in range(nr):
        elastic[:ts[r], r] = 0.0

    elastic.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/data/migration_test_elastic_iso_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")