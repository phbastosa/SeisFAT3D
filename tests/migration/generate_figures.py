import numpy as np
import matplotlib.pyplot as plt

from sys import path
path.append("../src/")
import functions

nx = 201
ny = 201
nz = 201

dx = 10.0
dy = 10.0
dz = 10.0

vp_model_file = "../inputs/models/diffraction_vp_201x201x201_10m.bin"   
vs_model_file = "../inputs/models/diffraction_vs_201x201x201_10m.bin"   
rho_model_file = "../inputs/models/diffraction_rho_201x201x201_10m.bin"   

shots_file = "../inputs/geometry/xyz_shots_position.txt"              
nodes_file = "../inputs/geometry/xyz_nodes_position.txt"  

shots = np.loadtxt(shots_file, dtype = np.float32, delimiter = ',')
nodes = np.loadtxt(nodes_file, dtype = np.float32, delimiter = ',')

model_vp = functions.read_binary_volume(nz, nx, ny, vp_model_file)
model_vs = functions.read_binary_volume(nz, nx, ny, vs_model_file)
model_rho = functions.read_binary_volume(nz, nx, ny, rho_model_file)

slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

dh = np.array([dx, dy, dz])

functions.plot_model_3D(model_vp, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 2000,
                        scale = 0.4, dbar = 1.25)

plt.savefig(f"vp_model_migration_test.png", dpi = 200)

functions.plot_model_3D(model_vs, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 2000,
                        scale = 0.4, dbar = 1.25)

plt.savefig(f"vs_model_migration_test.png", dpi = 200)

functions.plot_model_3D(model_rho, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 2000,
                        scale = 0.4, dbar = 1.25)

plt.savefig(f"rho_model_migration_test.png", dpi = 200)

# input data -------------------------------------------------------------------------------------

nt = 2001
nr = len(nodes)
ns = len(shots)

dt = 1e-3

fig, ax = plt.subplots(ncols = 1, nrows = 4, figsize = (8,10))

tloc = np.linspace(0, nt-1, 11)
tlab = np.around(tloc*dt, decimals = 1)

xloc = np.linspace(0, nr-1, 11, dtype = int)

for i in range(ns):
    
    input_data = functions.read_binary_matrix(nt, nr, f"../inputs/data/migration_test_data_nSamples{nt}_nRec{nr}_shot_{i+1}.bin")

    scale = np.std(input_data) 

    ax[i].imshow(input_data, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
    
    ax[i].set_xticks(xloc)
    ax[i].set_xticklabels(xloc)

    ax[i].set_yticks(tloc)
    ax[i].set_yticklabels(tlab)
    
    ax[i].set_title(f"Shot {i+1}", fontsize = 15)
    ax[i].set_xlabel(f"Trace index", fontsize = 15)
    ax[i].set_ylabel(f"Time [s]", fontsize = 15)


fig.tight_layout()
plt.savefig("input_seismograms_migration_test.png", dpi = 200)

# migrated data --------------------------------------------------------------------------------

image_file = "../outputs/images/image_kirchhoff_migration_test.bin"   

image = functions.read_binary_volume(nz, nx, ny, image_file)

functions.plot_model_3D(image, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        cmap = "seismic",
                        scale = 0.4, dbar = 1.25)

plt.savefig(f"image_migration_test.png", dpi = 200)


