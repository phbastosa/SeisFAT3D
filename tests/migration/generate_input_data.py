import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

parameters = str(sys.argv[1])

SPS = np.loadtxt(pyf.catch_parameter(parameters,"SPS"), delimiter = ",", dtype = np.float32) 
RPS = np.loadtxt(pyf.catch_parameter(parameters,"RPS"), delimiter = ",", dtype = np.float32) 
XPS = np.loadtxt(pyf.catch_parameter(parameters,"XPS"), delimiter = ",", dtype = np.int32) 

nt = int(pyf.catch_parameter(parameters, "time_samples"))
dt = float(pyf.catch_parameter(parameters, "time_spacing"))

fmax = 75.0

spread = XPS[0,2] - XPS[0,1]

vp = np.array([1500, 1700, 1900, 2300])
ro = np.array([1000, 2250, 2270, 2290])
z = np.array([300, 300, 300])

wavelet = pyf.get_ricker_wavelet(nt, dt, fmax)

input_data_folder = pyf.catch_parameter(parameters, "input_data_folder")
input_data_prefix = pyf.catch_parameter(parameters, "input_data_prefix")

for sId in range(len(SPS)):

    print(f"Creating data {sId+1} of {len(SPS)}")

    sx = SPS[sId,0]
    sz = SPS[sId,1]

    rx = RPS[XPS[sId,1]:XPS[sId,2],0]
    rz = RPS[XPS[sId,1]:XPS[sId,2],1]

    offset = np.sqrt((sx - rx)**2 + (rz - sz)**2)

    reflection_amps = pyf.get_analytical_amps_reflections(vp, ro, z)
    reflection_time = pyf.get_analytical_time_reflections(vp, z, offset)

    input_data = np.zeros((nt, spread))

    for trace in range(spread):
        for layer in range(len(z)):
            time_index = int(reflection_time[layer,trace] / dt)

            if time_index < nt:
                input_data[time_index,trace] = reflection_amps[layer]

        input_data[:,trace] = np.convolve(input_data[:,trace], wavelet, "same") 

    input_data.flatten("F").astype(np.float32, order = "F").tofile(input_data_folder + input_data_prefix + f"{sId+1}.bin")
    