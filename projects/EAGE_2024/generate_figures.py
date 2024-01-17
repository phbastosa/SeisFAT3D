import numpy as np
import matplotlib.pyplot as plt

from sys import path
path.append("../src/")
import functions

#---------------------------------------------------------------------

nx = 441
ny = 441
nz = 161

dh = np.array([25, 25, 25])

model = functions.read_binary_volume(nz, nx, ny, f"../inputs/models/accuracyModelTest_{nz}x{nx}x{ny}_{dh[0]:.0f}m.bin")

shots_file = "../inputs/geometry/xyz_shots_position.txt"
nodes_file = "../inputs/geometry/xyz_nodes_position.txt"

slices = np.array([int(nz/2), int(nx/2), 240], dtype = int) # [xy, zy, zx]

eikonal = np.zeros((3, nz, nx, ny))

eikonal[0] = functions.read_binary_volume(nz, nx, ny, f"../outputs/travel_times/fim_time_volume_{nz}x{nx}x{ny}_shot_1.bin")
eikonal[1] = functions.read_binary_volume(nz, nx, ny, f"../outputs/travel_times/fsm_time_volume_{nz}x{nx}x{ny}_shot_1.bin")
eikonal[2] = functions.read_binary_volume(nz, nx, ny, f"../outputs/travel_times/hd_fim_time_volume_{nz}x{nx}x{ny}_shot_1.bin")

eikonal_levels = [1.0, 2.0, 3.0]
eikonal_colors = ["orange", "blue", "green"]

functions.plot_model_3D(model, dh, slices,
                        shots = shots_file, 
                        nodes = nodes_file, 
                        eikonal = eikonal,
                        eikonal_colors = eikonal_colors,
                        eikonal_levels = eikonal_levels,
                        scale = 1.6)

plt.savefig(f"accuracy_model.png", dpi = 200)
plt.clf()

#---------------------------------------------------------------------

shots = np.loadtxt(shots_file, delimiter = ',') 
nodes = np.loadtxt(nodes_file, delimiter = ',') 

nTraces = len(nodes)

velocity = np.array([2000, 3000, 4000, 5000])
thickness = np.array([1000, 1000, 1000])

offset = np.sqrt((shots[0] - nodes[:,0])**2 + (shots[1] - nodes[:,1])**2 + (shots[2] - nodes[:,2])**2)

analytical_times = functions.analytical_first_arrivals(velocity, thickness, offset)

fsm_firstArrivals = np.fromfile(f"../outputs/first_arrivals/fsm_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)
fim_firstArrivals = np.fromfile(f"../outputs/first_arrivals/fim_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)
ifim_firstArrivals = np.fromfile(f"../outputs/first_arrivals/hd_fim_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)

sl1 = slice(0, int(nTraces/3))
sl2 = slice(int(nTraces/3), int(2*nTraces/3))
sl3 = slice(int(2*nTraces/3), nTraces)

slices = [sl1, sl2, sl3]

limits = [[2, 3.5], [2.5, 4.0], [3.0, 4.5]]

xloc = np.linspace(0, nTraces/3-1, 7)
xlab = np.linspace(0, nTraces/3, 7, dtype = int)

fig, ax = plt.subplots(ncols = 3, nrows = 1, figsize = (10, 4))

for i in range(len(slices)):

    tloc = np.linspace(limits[i][0], limits[i][1], 5)
    tlab = np.around(tloc, decimals = 1)

    ax[i].plot(analytical_times[slices[i]], "-", color = "red", markersize = 1, label = "Analytical travel times")
    ax[i].plot(fsm_firstArrivals[slices[i]], "--", color = "blue", markersize = 1, label = "Noble, Gesret & Belayouni (2014)")
    ax[i].plot(fim_firstArrivals[slices[i]], "--", color = "orange", markersize = 1, label = "Jeong & Whitaker (2008)")
    ax[i].plot(ifim_firstArrivals[slices[i]], "--", color = "green", markersize = 1, label = "Cai, Zhu & Li (2023)")

    ax[i].set_yticks(tloc)
    ax[i].set_yticklabels(tlab)

    ax[i].set_xticks(xloc)
    ax[i].set_xticklabels(xlab)

    ax[i].set_ylim(limits[i])

    ax[i].invert_yaxis()

    ax[i].set_xlabel("Trace index", fontsize = 12)
    ax[i].set_ylabel("Time [s]", fontsize = 12)

    ax[i].legend(loc = "upper right", fontsize = 8)

fig.tight_layout()
plt.savefig("seismograms.png", dpi = 300)

limits = [[-0.15, 0.1], [-0.15, 0.1], [-0.15, 0.1]]

fig, ax = plt.subplots(ncols = 3, nrows = 1, figsize = (10, 4))

for i in range(len(slices)):

    tloc = np.linspace(limits[i][0], limits[i][1], 11)
    tlab = np.around(tloc, decimals = 2)

    ax[i].plot(analytical_times[slices[i]] - fsm_firstArrivals[slices[i]], "--", color = "blue", markersize = 1, label = "Noble, Gesret & Belayouni (2014)")
    ax[i].plot(analytical_times[slices[i]] - fim_firstArrivals[slices[i]], "--", color = "orange", markersize = 1, label = "Jeong & Whitaker (2008)")
    ax[i].plot(analytical_times[slices[i]] - ifim_firstArrivals[slices[i]], "--", color = "green", markersize = 1, label = "Cai, Zhu & Li (2023)")

    ax[i].set_yticks(tloc)
    ax[i].set_yticklabels(tlab)

    ax[i].set_xticks(xloc)
    ax[i].set_xticklabels(xlab)

    ax[i].set_ylim(limits[i])

    ax[i].set_xlabel("Trace index", fontsize = 12)
    ax[i].set_ylabel("Time difference [s]", fontsize = 12)

    ax[i].legend(loc = "upper right", fontsize = 8)

fig.tight_layout()
plt.savefig("diff_seismograms.png", dpi = 300)

max_error_fsm = np.max(analytical_times - fsm_firstArrivals)
min_error_fsm = np.min(analytical_times - fsm_firstArrivals)
std_error_fsm = np.std(analytical_times - fsm_firstArrivals)
mean_error_fsm = np.mean(analytical_times - fsm_firstArrivals)

max_error_fim = np.max(analytical_times - fim_firstArrivals)
min_error_fim = np.min(analytical_times - fim_firstArrivals)
std_error_fim = np.std(analytical_times - fim_firstArrivals)
mean_error_fim = np.mean(analytical_times - fim_firstArrivals)

max_error_ifim = np.max(analytical_times - ifim_firstArrivals)
min_error_ifim = np.min(analytical_times - ifim_firstArrivals)
std_error_ifim = np.std(analytical_times - ifim_firstArrivals)
mean_error_ifim = np.mean(analytical_times - ifim_firstArrivals)

print(f"|{'-'*98}|")
print(f"| {'Data difference (ta - tn)':^40s} | {'MAX ERROR':^11s} | {'MIN ERROR':^11s} | {'STD ERROR':^11s} | {'MEAN ERROR':^11s} |")
print(f"|{'-'*98}|")
print(f"| {'Noble, Gesret & Belayouni (2014)':^40s} | {f'{max_error_fsm:.4f}':^11s} | {f'{min_error_fsm:.5f}':^11s} | {f'{std_error_fsm:.4f}':^11s} | {f'{mean_error_fsm:.4f}':^11s} |")
print(f"| {'Jeong & Whitaker (2008)':^40s} | {f'{max_error_fim:.4f}':^11s} | {f'{min_error_fim:.5f}':^11s} | {f'{std_error_fim:.4f}':^11s} | {f'{mean_error_fim:.4f}':^11s} |")
print(f"| {'Cai, Zhu & Li (2023)':^40s} | {f'{max_error_ifim:.4f}':^11s} | {f'{min_error_ifim:.5f}':^11s} | {f'{std_error_ifim:.4f}':^11s} | {f'{mean_error_ifim:.4f}':^11s} |")
print(f"|{'-'*98}|\n\n")
