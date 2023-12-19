import numpy as np
import matplotlib.pyplot as plt

from sys import path

path.append("../src/")

import functions

#---------------------------------------------------------------------

nx = 441
ny = 441
nz = 161

dh = 25

accuracy_model = functions.read_binary_volume(nz, nx, ny, f"../inputs/models/accuracyModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

shots_file = "../inputs/geometry/xyz_shots_position.txt"
nodes_file = "../inputs/geometry/xyz_nodes_position.txt"
 
shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([int(nz/2), int(nx/2), 240], dtype = int) # [xy, zy, zx]
dh = np.array([dh, dh, dh])

vmin = 2000
vmax = 5000

volFIM = functions.read_binary_volume(nz, nx, ny, f"../outputs/travel_times/fim_time_volume_{nz}x{nx}x{ny}_shot_1.bin")
volFSM = functions.read_binary_volume(nz, nx, ny, f"../outputs/travel_times/fsm_time_volume_{nz}x{nx}x{ny}_shot_1.bin")

functions.plot_model_eikonal_3D(accuracy_model, volFIM, volFSM, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
plt.savefig(f"accuracy_model_test.png", dpi = 200)
plt.clf()

nt = 6001
dt = 1e-3

nTraces = len(nodes)

velocity = np.array([2000, 3000, 4000, 5000])
thickness = np.array([1000, 1000, 1000])

offset = np.sqrt((shots[0] - nodes[:,0])**2 + (shots[1] - nodes[:,1])**2 + (shots[2] - nodes[:,2])**2)

analytical_times = functions.analytical_first_arrivals(velocity, thickness, offset)

fim_firstArrivals = np.fromfile(f"../outputs/first_arrivals/fim_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)
fsm_firstArrivals = np.fromfile(f"../outputs/first_arrivals/fsm_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)

t0 = int(np.pi / 25.0 / dt) 
tmax = int(5.0 / dt) + 1

seismic = functions.read_binary_matrix(nt, nTraces, f"../inputs/data/synthetic_seismogram_6001x273_shot_1.bin")

seismic = seismic[t0:tmax+t0,:]

sl1 = slice(0, int(nTraces/3))
sl2 = slice(int(nTraces/3), int(2*nTraces/3))
sl3 = slice(int(2*nTraces/3), nTraces)

slices = [sl1, sl2, sl3]

scale = 0.1*np.std(seismic)

tloc = np.linspace(0, tmax-1, 11)
tlab = np.around(tloc * dt, decimals = 3)

xloc = np.linspace(0, nTraces/3-1, 7)
xlab = np.linspace(0, nTraces/3, 7, dtype = int)

fig, ax = plt.subplots(ncols = 3, nrows = 1, figsize = (10,6))

for i in range(len(slices)):

    ax[i].imshow(seismic[:, slices[i]], aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
    ax[i].plot(analytical_times[slices[i]]/dt, "o", color = "red", markersize = 1, label = "Analytical travel times")
    ax[i].plot(fim_firstArrivals[slices[i]]/dt, "o", color = "orange", markersize = 1, label = "Jeong & Whitaker (2008)")
    ax[i].plot(fsm_firstArrivals[slices[i]]/dt, "o", color = "green", markersize = 1, label = "Noble, Gesret & Belayouni (2014)")

    ax[i].set_yticks(tloc)
    ax[i].set_yticklabels(tlab)

    ax[i].set_xticks(xloc)
    ax[i].set_xticklabels(xlab)

    ax[i].set_xlabel("Trace index", fontsize = 12)
    ax[i].set_ylabel("Time [s]", fontsize = 12)

    ax[i].legend(loc = "upper right", fontsize = 8)

fig.tight_layout()
plt.savefig("modeled_seismograms.png", dpi = 300)

