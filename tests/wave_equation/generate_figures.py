import numpy as np
import matplotlib.pyplot as plt

from os import system

from sys import path, argv
path.append("../src/")
import functions

nx = int(functions.catch_parameter(argv[1], "x_samples"))
ny = int(functions.catch_parameter(argv[1], "y_samples"))
nz = int(functions.catch_parameter(argv[1], "z_samples"))

dx = float(functions.catch_parameter(argv[1], "x_spacing"))
dy = float(functions.catch_parameter(argv[1], "y_spacing"))
dz = float(functions.catch_parameter(argv[1], "z_spacing"))

vp_model_file = functions.catch_parameter(argv[1], "vp_model_file")
vs_model_file = functions.catch_parameter(argv[1], "vs_model_file")
rho_model_file = functions.catch_parameter(argv[1], "rho_model_file")

shots_file = functions.catch_parameter(argv[1], "shots_file")
nodes_file = functions.catch_parameter(argv[1], "nodes_file")

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
                        scale = 0.35, dbar = 1.3)

plt.savefig(f"vp_model_test.png", dpi = 200)

functions.plot_model_3D(model_vs, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 3000,
                        scale = 0.35, dbar = 1.3)

plt.savefig(f"vs_model_test.png", dpi = 200)

functions.plot_model_3D(model_rho, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 3000,
                        scale = 0.35, dbar = 1.3)

plt.savefig(f"rho_model_test.png", dpi = 200)

nTraces = len(nodes)

nt = int(functions.catch_parameter(argv[1], "time_samples"))
dt = float(functions.catch_parameter(argv[1], "time_spacing"))

fmax = float(functions.catch_parameter(argv[1], "max_frequency"))

tlag = 2.0*np.sqrt(np.pi) / fmax

scalar_seismogram_path = f"../outputs/seismograms/scalar_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"
elastic_seismogram_path = f"../outputs/seismograms/elastic_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"
acoustic_seismogram_path = f"../outputs/seismograms/acoustic_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"

scalar_seism = functions.read_binary_matrix(nt, nTraces, scalar_seismogram_path)
elastic_seism = functions.read_binary_matrix(nt, nTraces, elastic_seismogram_path)
acoustic_seism = functions.read_binary_matrix(nt, nTraces, acoustic_seismogram_path)

trace_index = int(0.5*nTraces)

scalar_time = dt*np.where(scalar_seism[:,trace_index] == np.max(scalar_seism[:,trace_index]))[0][0]
acoustic_time = dt*np.where(acoustic_seism[:,trace_index] == np.max(acoustic_seism[:,trace_index]))[0][0]
elastic_time = dt*np.where(elastic_seism[:,trace_index] == np.max(elastic_seism[:,trace_index]))[0][0]

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
ax[3].set_title("Refractions", fontsize = 15)
ax[3].set_xlabel("Amplitude", fontsize = 15)
ax[3].set_ylabel("Time [s]", fontsize = 15)
ax[3].set_ylim([tmin, tmax])
ax[3].legend(loc = "upper right")
ax[3].invert_yaxis()

fig.tight_layout()
plt.savefig("seismograms.png", dpi = 200)