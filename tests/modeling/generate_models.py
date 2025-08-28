import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

parameters = str(sys.argv[1])

nx = int(pyf.catch_parameter(parameters, "x_samples"))
ny = int(pyf.catch_parameter(parameters, "y_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples"))

dz = float(pyf.catch_parameter(parameters, "z_spacing"))

vp = np.array([1500, 1700, 1900, 2300, 3000, 3500])
z = np.array([200, 500, 1000, 1500, 1500])

Vp = np.zeros((nz,nx,ny)) + vp[0]
for i in range(len(z)): 
    Vp[int(np.sum(z[:i+1]/dz)):] = vp[i+1]

Vp.flatten("F").astype(np.float32, order = "F").tofile(pyf.catch_parameter(parameters, "vp_model_file"))
