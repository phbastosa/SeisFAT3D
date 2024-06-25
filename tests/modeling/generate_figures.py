import numpy as np
import matplotlib.pyplot as plt

from sys import path, argv
path.append("../src/")
import functions

#-------------------------------------------------------------------------

nx = int(functions.catch_parameter(argv[1], "x_samples"))
ny = int(functions.catch_parameter(argv[1], "y_samples"))
nz = int(functions.catch_parameter(argv[1], "z_samples"))

dx = float(functions.catch_parameter(argv[1], "x_spacing"))
dy = float(functions.catch_parameter(argv[1], "y_spacing"))
dz = float(functions.catch_parameter(argv[1], "z_spacing"))

model = functions.read_binary_volume(nz, nx, ny, functions.catch_parameter(argv[1], "vp_model_file"))

shots_file = functions.catch_parameter(argv[1], "shots_file")
nodes_file = functions.catch_parameter(argv[1], "nodes_file")

slices = np.array([nz/2, ny/2, nx/2], dtype = int)
dh = np.array([dx, dy, dz])

functions.plot_model_3D(model, dh, slices, 
                        shots = shots_file, 
                        nodes = nodes_file,
                        scale = 2.5, 
                        dbar = 1.6)

plt.savefig(f"test_model.png", dpi = 300)

#-------------------------------------------------------------------------

nt = int(functions.catch_parameter(argv[1], "time_samples"))
dt = float(functions.catch_parameter(argv[1], "time_spacing"))

fmax = float(functions.catch_parameter(argv[1], "max_frequency"))

tlag = 2.0*np.sqrt(np.pi)/fmax

shots = np.loadtxt(shots_file, delimiter = ",", dtype = float)
nodes = np.loadtxt(nodes_file, delimiter = ",", dtype = float)

nTraces = len(nodes)

eikonal = functions.read_binary_array(nTraces, f"../outputs/seismograms/eikonal_data_nRec{nTraces}_shot_1.bin")
seismic = functions.read_binary_matrix(nt, nTraces, f"../outputs/seismograms/scalar_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin")

scale = 0.01*np.std(seismic)

tloc = np.linspace(0, nt, 11)
tlab = np.around(tloc * dt, decimals = 1)

fig, ax = plt.subplots(figsize = (15,5))

ax.imshow(seismic[:,:500], aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
ax.plot(eikonal[:500]/dt + tlag, color = "red")

ax.set_yticks(tloc)
ax.set_yticklabels(tlab)

plt.tight_layout()
plt.savefig("modeling_results.png", dpi = 300)