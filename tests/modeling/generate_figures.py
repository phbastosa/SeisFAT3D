import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

parameters = str(sys.argv[1])

sps_path = pyf.catch_parameter(parameters, "SPS") 
rps_path = pyf.catch_parameter(parameters, "RPS") 
xps_path = pyf.catch_parameter(parameters, "XPS") 

nx = int(pyf.catch_parameter(parameters, "x_samples"))
ny = int(pyf.catch_parameter(parameters, "y_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples"))

dx = float(pyf.catch_parameter(parameters, "x_spacing"))
dy = float(pyf.catch_parameter(parameters, "y_spacing"))
dz = float(pyf.catch_parameter(parameters, "z_spacing"))

model_vp = pyf.read_binary_volume(nz, nx, ny, pyf.catch_parameter(parameters, "vp_model_file"))

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model_vp, dh, slices, shots = sps_path, nodes = rps_path, scale = 14, 
                 adjx = 0.8, dbar = 1.5, cmap = "Greys", cblab = "Vp [m/s]")
plt.savefig("modeling_test_setup.png", dpi = 200)
plt.show()

SPS = np.loadtxt(sps_path, delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt(rps_path, delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt(xps_path, delimiter = ",", comments = "#", dtype = int)  

ns = len(SPS)
nr = len(RPS)

v = np.array([1500, 1700, 1900, 2300, 3000, 3500])
z = np.array([200, 500, 1000, 1500, 1500])

output_folder = pyf.catch_parameter(parameters, "modeling_output_folder")

fig, ax = plt.subplots(figsize = (15, 7), nrows = 3)

eikonal_an = np.zeros(nr)

for i in range(ns):

    x = np.sqrt((SPS[i,0] - RPS[:,0])**2 + (SPS[i,1] - RPS[:,1])**2)

    refractions = pyf.get_analytical_refractions(v,z,x)

    for k in range(nr):
        eikonal_an[k] = min(x[k]/v[0], np.min(refractions[:,k]))
    
    eikonal_nu = pyf.read_binary_array(nr, output_folder + f"eikonal_iso_nStations{nr}_shot_{i+1}.bin")

    ax[i].plot(eikonal_an - eikonal_nu, "k")

    ax[i].set_ylabel("(Ta - Tn) [ms]", fontsize = 15)
    ax[i].set_xlabel("Channel index", fontsize = 15)
    
    ax[i].set_yticks(np.linspace(-0.005, 0.005, 5))
    ax[i].set_yticklabels(np.linspace(-5, 5, 5, dtype = float))

    ax[i].set_xticks(np.linspace(0, nr, 11))
    ax[i].set_xticklabels(np.linspace(0, nr, 11, dtype = int))

    ax[i].set_xlim([0, nr])

    ax[i].invert_yaxis()

fig.tight_layout()
plt.savefig("modeling_test_accuracy.png", dpi = 200)
plt.show()