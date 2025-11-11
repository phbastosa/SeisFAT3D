import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

parameters = str(sys.argv[1])

SPS_path = pyf.catch_parameter(parameters, "SPS")
RPS_path = pyf.catch_parameter(parameters, "RPS")
XPS_path = pyf.catch_parameter(parameters, "XPS")

SPS = np.loadtxt(SPS_path, delimiter = ",", dtype = np.float32) 
RPS = np.loadtxt(RPS_path, delimiter = ",", dtype = np.float32) 
XPS = np.loadtxt(XPS_path, delimiter = ",", dtype = np.int32) 

nt = int(pyf.catch_parameter(parameters, "time_samples"))
dt = float(pyf.catch_parameter(parameters, "time_spacing"))

ns = len(SPS)
nr = len(RPS)

spread = XPS[0,2]

for sId in range(ns):

    data_file = f"../../FWI3D/outputs/data/seismogram_nt4001_nr1200_1000us_shot_{sId+1}.bin"

    data = pyf.read_binary_matrix(nt, spread, data_file)

    data *= 1.0 / np.max(np.abs(data))

    travel_time = np.sqrt((SPS[sId,0] - RPS[XPS[sId,1]:XPS[sId,2],0])**2 + 
                          (SPS[sId,1] - RPS[XPS[sId,1]:XPS[sId,2],1])**2) / 1500

    tId = np.array((travel_time + 0.15) / dt, dtype = int)

    sigma = 50
    for rId in range(spread):
        if tId[rId] >= nt:
            tId[rId] = nt-1 

        data[:tId[rId], rId] *= np.exp(-((np.arange(tId[rId])-tId[rId]) / sigma)**2)        

    data[:900,:] = 0.0

    output_path = f"../inputs/data/seismogram_input_shot_{sId+1}.bin"
    
    data.flatten("F").astype(np.float32, order = "F").tofile(output_path)

import matplotlib.pyplot as plt

plt.imshow(data, aspect = "auto", cmap = "Greys")
plt.show()
