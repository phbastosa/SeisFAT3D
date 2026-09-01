import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf
import matplotlib.pyplot as plt

parameters = str(sys.argv[1])

sps_path = pyf.catch_parameter(parameters,"SPS") 
rps_path = pyf.catch_parameter(parameters,"RPS") 

nx = int(pyf.catch_parameter(parameters, "x_samples"))
ny = int(pyf.catch_parameter(parameters, "y_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples")) 

dh = float(pyf.catch_parameter(parameters, "model_spacing"))

image = np.fromfile(f"../outputs/seismic/IDKDM_result_161x401x401.bin", count = nx*ny*nz, dtype = np.float32).reshape([nz,nx,ny], order = "F")

image *= -2000 / np.max(np.abs(image))

SPS = np.loadtxt(sps_path, delimiter = ",", dtype = float)
RPS = np.loadtxt(rps_path, delimiter = ",", dtype = float)

slices = np.array([700/dh, 0.5*ny, 1000/dh], dtype = int)

dh = np.array([dh, dh, dh])

pyf.plot_model_3D(image, dh, slices, shots = sps_path, 
                  scale = 1.4, adjx = 0.7, dbar = 1.4, cmap = "Greys",
                  cblab = "Normalized Amplitude", vmin = -500, vmax = 500)
plt.show()

