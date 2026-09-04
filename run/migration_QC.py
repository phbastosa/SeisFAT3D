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

IDKDM1 = np.fromfile(f"first_src_line.bin", count = nx*ny*nz, dtype = np.float32).reshape([nz,nx,ny], order = "F")
IDKDM2 = np.fromfile(f"last_src_line.bin", count = nx*ny*nz, dtype = np.float32).reshape([nz,nx,ny], order = "F")

IDRTM = np.fromfile(f"../../FWI3D/outputs/seismic/RTM_section_30Hz_161x401x401_12m.bin", count = nx*ny*nz, dtype = np.float32).reshape([nz,nx,ny], order = "F")

IDKDM1 *= 2000 / np.max(np.abs(IDKDM1))
IDKDM2 *= 2000 / np.max(np.abs(IDKDM2))

IDKDM = IDKDM1 + IDKDM2

IDRTM *= 2000 / np.max(np.abs(IDRTM))

SPS = np.loadtxt(sps_path, delimiter = ",", dtype = float)
RPS = np.loadtxt(rps_path, delimiter = ",", dtype = float)

slices = np.array([700/dh, 0.5*ny, 0.5*nx], dtype = int)

dh = np.array([dh, dh, dh])

pyf.plot_model_3D(IDKDM, dh, slices, shots = sps_path, 
                  scale = 1.4, adjx = 0.7, dbar = 1.4, cmap = "Greys",
                  cblab = "Normalized Amplitude", vmin = -500, vmax = 500)
#plt.show()

pyf.plot_model_3D(IDRTM, dh, slices, shots = sps_path, 
                  scale = 1.4, adjx = 0.7, dbar = 1.4, cmap = "Greys",
                  cblab = "Normalized Amplitude", vmin = -500, vmax = 500)
plt.show()
