import sys.path; path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

from os import system

nx = 201
ny = 201
nz = 201

dx = 10.0
dy = 10.0
dz = 10.0

vp_model_file = "../inputs/models/homogeneous_vp_201x201x201_10m.bin"   
vs_model_file = "../inputs/models/homogeneous_vs_201x201x201_10m.bin"   
rho_model_file = "../inputs/models/homogeneous_rho_201x201x201_10m.bin"   

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
                        vmax = 3000,
                        scale = 0.4, dbar = 1.25)

plt.savefig(f"vp_model_wave_equation_test.png", dpi = 200)

functions.plot_model_3D(model_vs, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 3000,
                        scale = 0.4, dbar = 1.25)

plt.savefig(f"vs_model_wave_equation_test.png", dpi = 200)

functions.plot_model_3D(model_rho, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 3000,
                        scale = 0.4, dbar = 1.25)

plt.savefig(f"rho_model_wave_equation_test.png", dpi = 200)

# Snapshots ------------------------------------------------------------------------------------

scalar_snapshots_path = "../outputs/snapshots/scalar_snapshot_201x201x201_shot_1_Nsnaps10.bin"
acoustic_snapshots_path = "../outputs/snapshots/acoustic_snapshot_201x201x201_shot_1_Nsnaps10.bin"
elastic_snapshots_path = "../outputs/snapshots/elastic_snapshot_201x201x201_shot_1_Nsnaps10.bin"

scalar_snaps = functions.read_binary_volume(nz, nx, 10*ny, scalar_snapshots_path)
acoustic_snaps = functions.read_binary_volume(nz, nx, 10*ny, acoustic_snapshots_path)
elastic_snaps = functions.read_binary_volume(nz, nx, 10*ny, elastic_snapshots_path)

functions.plot_model_3D(scalar_snaps[:,:,3*ny:3*ny+ny], dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        cmap = "Greys",
                        scale = 0.4, 
                        dbar = 1.25)

plt.savefig(f"scalar_snapshot_wave_equation_test.png", dpi = 200)

functions.plot_model_3D(acoustic_snaps[:,:,3*ny:3*ny+ny], dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        cmap = "Greys",
                        scale = 0.4, 
                        dbar = 1.25)

plt.savefig(f"acoustic_snapshot_wave_equation_test.png", dpi = 200)

functions.plot_model_3D(elastic_snaps[:,:,3*ny:3*ny+ny], dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        cmap = "Greys",
                        scale = 0.4, 
                        dbar = 1.25)

plt.savefig(f"elastic_snapshot_wave_equation_test.png", dpi = 200)

# Seismograms ----------------------------------------------------------------------------------

nTraces = len(nodes)

nt = 3001
dt = 1e-3

fmax = 30.0

tlag = 2.0*np.sqrt(np.pi) / fmax

scalar_seismogram_path = f"../outputs/seismograms/scalar_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"
elastic_seismogram_path = f"../outputs/seismograms/elastic_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"
acoustic_seismogram_path = f"../outputs/seismograms/acoustic_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"

scalar_seism = functions.read_binary_matrix(nt, nTraces, scalar_seismogram_path)
elastic_seism = functions.read_binary_matrix(nt, nTraces, elastic_seismogram_path)
acoustic_seism = functions.read_binary_matrix(nt, nTraces, acoustic_seismogram_path)

trace_index = int(0.5*nTraces)

scalar_time = dt*(np.where(scalar_seism[:,trace_index] == np.max(scalar_seism[:,trace_index]))[0][0]+1)
acoustic_time = dt*(np.where(acoustic_seism[:,trace_index] == np.max(acoustic_seism[:,trace_index]))[0][0]+1)
elastic_time = dt*(np.where(elastic_seism[:,trace_index] == np.max(elastic_seism[:,trace_index]))[0][0]+1)

system("clear")
print(f"Analytical travel time = {1800/1500 + tlag} s")
print(f"Scalar travel time = {scalar_time} s")
print(f"Acoustic travel time = {acoustic_time} s")
print(f"Elastic travel time = {elastic_time} s")

scale = 10.0*np.std(elastic_seism)

tloc = np.linspace(0, nt-1, 11)
tlab = np.around(tloc * dt, decimals = 1)

xloc = np.linspace(0, nTraces, 5, dtype = int)

tmin = 1.1
tmax = 1.5
xmin = 0.5*nTraces - 5
xmax = 0.5*nTraces + 5

rectangle = np.array([[tmax, xmin], [tmin, xmin], [tmin, xmax], [tmax, xmax], [tmax, xmin]])

fig, ax = plt.subplots(ncols = 4, nrows = 1, figsize = (15,6))

im1 = ax[0].imshow(scalar_seism, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
ax[0].set_title("Scalar media", fontsize = 15)
cbar1 = fig.colorbar(im1, ax = ax[0])
cbar1.set_label("amplitude", fontsize = 12)
ax[0].set_yticks(tloc)
ax[0].set_yticklabels(tlab)
ax[0].set_xticks(xloc)
ax[0].set_xticklabels(xloc)
ax[0].set_xlabel("Trace number", fontsize = 15)
ax[0].set_ylabel("Time [s]", fontsize = 15)
ax[0].plot(rectangle[:,1], rectangle[:,0]/dt)

im2 = ax[1].imshow(acoustic_seism, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
ax[1].set_title("Acoustic media", fontsize = 15)
cbar2 = fig.colorbar(im2, ax = ax[1])
cbar2.set_label("amplitude", fontsize = 12)
ax[1].set_yticks(tloc)
ax[1].set_yticklabels(tlab)
ax[1].set_xticks(xloc)
ax[1].set_xticklabels(xloc)
ax[1].set_xlabel("Trace number", fontsize = 15)
ax[1].set_ylabel("Time [s]", fontsize = 15)
ax[1].plot(rectangle[:,1], rectangle[:,0]/dt)

im3 = ax[2].imshow(elastic_seism, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
ax[2].set_title("Elastic media", fontsize = 15)
cbar3 = fig.colorbar(im3, ax = ax[2])
cbar3.set_label("amplitude", fontsize = 12)
ax[2].set_yticks(tloc)
ax[2].set_yticklabels(tlab)
ax[2].set_xticks(xloc)
ax[2].set_xticklabels(xloc)
ax[2].set_xlabel("Trace number", fontsize = 15)
ax[2].set_ylabel("Time [s]", fontsize = 15)
ax[2].plot(rectangle[:,1], rectangle[:,0]/dt)

time = np.arange(nt)*dt
mask = np.logical_and(time > tmin, time < tmax)

ax[3].plot(scalar_seism[mask,int(0.5*nTraces)], time[mask], label = "Scalar media")
ax[3].plot(acoustic_seism[mask,int(0.5*nTraces)], time[mask], label = "Acoustic media")
ax[3].plot(elastic_seism[mask,int(0.5*nTraces)], time[mask], label = "Elastic media")
ax[3].set_title("Central trace", fontsize = 15)
ax[3].set_xlabel("Amplitude", fontsize = 15)
ax[3].set_ylabel("Time [s]", fontsize = 15)
ax[3].set_ylim([tmin, tmax])
ax[3].legend(loc = "upper right")
ax[3].invert_yaxis()

fig.tight_layout()
plt.savefig("seismograms_from_wave_equation_test.png", dpi = 200)