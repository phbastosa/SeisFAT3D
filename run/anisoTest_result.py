import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

path_SPS = "../inputs/geometry/anisoTest_SPS.txt"
path_RPS = "../inputs/geometry/anisoTest_RPS.txt"
path_XPS = "../inputs/geometry/anisoTest_XPS.txt"

SPS = np.loadtxt(path_SPS, dtype = float, delimiter = ",")
RPS = np.loadtxt(path_RPS, dtype = float, delimiter = ",")
XPS = np.loadtxt(path_XPS, dtype = float, delimiter = ",")

nx = 201
ny = 201
nz = 101

dx = 10.0
dy = 10.0
dz = 10.0

model = pyf.read_binary_volume(nz, nx, ny, "../inputs/models/anisoTest_vp.bin")

dh = np.array([dx, dy, dz])
slices = np.array([0.5*nz, 0.5*ny, 0.5*nx], dtype = int)

pyf.plot_model_3D(model, dh, slices, shots = path_SPS, 
                  nodes = path_RPS, scale = 1.1, adjx = 0.65, 
                  dbar = 1.3, cblab = "P wave velocity [km/s]")
plt.show()

nt = 2001
dt = 1e-3
nr = 164

elastic_iso = pyf.read_binary_matrix(nt,nr,"../outputs/syntheticData/elastic_iso_nStations164_nSamples2001_shot_1.bin")
elastic_ani = pyf.read_binary_matrix(nt,nr,"../outputs/syntheticData/elastic_ani_nStations164_nSamples2001_shot_1.bin")

scale = 5.0*np.std(elastic_iso)

fig, ax = plt.subplots(nrows = 2, figsize = (12,8))

ax[0].imshow(elastic_iso, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale, extent = [0, 164, 2, 0])
ax[0].set_xlabel("Receiver index", fontsize = 15)
ax[0].set_ylabel("Time [s]", fontsize = 15)
ax[0].set_title("SSG", fontsize = 15)

ax[1].imshow(elastic_ani, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale, extent = [0, 164, 2, 0])
ax[1].set_xlabel("Receiver index", fontsize = 15)
ax[1].set_ylabel("Time [s]", fontsize = 15)
ax[1].set_title("RSG", fontsize = 15)

plt.tight_layout()
plt.savefig("gather.png", dpi = 300)
plt.show()

trace_iso = elastic_iso[:,60]
trace_ani = elastic_ani[:,60]

time = np.arange(nt)*dt

freq = np.fft.fftfreq(nt, dt)
signal_iso = np.fft.fft(trace_iso)
signal_ani = np.fft.fft(trace_ani)

mask = freq >= 0 

fig, ax = plt.subplots(ncols=2, figsize = (6,8))

ax[0].plot(trace_iso, time, label = "SSG")
ax[0].plot(trace_ani, time, label = "RSG")
ax[0].set_ylim([0, 2.0])
ax[0].invert_yaxis()
ax[0].set_xlabel("Amplitude", fontsize = 15)
ax[0].set_ylabel("Time [s]", fontsize = 15)
ax[0].legend(loc = "lower right")

ax[1].plot(np.abs(signal_iso[mask]), freq[mask], label = "SSG")
ax[1].plot(np.abs(signal_ani[mask]), freq[mask], label = "RSG")
ax[1].set_ylim([0, 60])
ax[1].invert_yaxis()
ax[1].set_xlabel("Amplitude", fontsize = 15)
ax[1].set_ylabel("Frequency [Hz]", fontsize = 15)
ax[1].legend(loc = "lower right")

plt.tight_layout()
plt.savefig("trace.png", dpi = 300)
plt.show()

# vp = np.array([2000, 3000])
# vs = np.array([1180, 1460])
# ro = np.array([2250, 2350])
# z = np.array([800])

# offset = np.sqrt((SPS[0] - RPS[:,0])**2 + (SPS[1] - RPS[:,1])**2)

# directWave = offset / vp[0] 
# refracions = pyf.get_analytical_refractions(vp, z, offset)[0]

# eikonal_iso = pyf.read_binary_array(nrx*nry,"../outputs/syntheticData/eikonal_iso_nStations603_shot_1.bin")
# eikonal_ani = pyf.read_binary_array(nrx*nry,"../outputs/syntheticData/eikonal_ani_nStations603_shot_1.bin")

# fig, ax = plt.subplots(ncols = 3, figsize = (15,9))

# for n in range(nrx):
#     ax[n].plot(directWave[n*nry:n*nry+nry], "o", color = "blue")
#     ax[n].plot(refracions[n*nry:n*nry+nry], "o", color = "black")
    
#     ax[n].plot(eikonal_iso[n*nry:n*nry+nry], "o", color = "green")
#     ax[n].plot(eikonal_ani[n*nry:n*nry+nry], "o", color = "yellow")

#     ax[n].set_ylim([0,10])
#     ax[n].invert_yaxis()
    
# plt.tight_layout()
# plt.show()
# sys.exit()




