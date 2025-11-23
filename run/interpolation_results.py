import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

nb = 3

nx = 801
ny = 801
nz = 171

dh1 = 25.0
dh2 = 50.0
dh3 = 100.0

nrx = 721
nry = 721

nx1 = int((nx-1)*dh1/dh1+1)
ny1 = int((ny-1)*dh1/dh1+1)
nz1 = int((nz-1)*dh1/dh1+1)

nx2 = int((nx-1)*dh1/dh2+1)
ny2 = int((ny-1)*dh1/dh2+1)
nz2 = int((nz-1)*dh1/dh2+1)

nx3 = int((nx-1)*dh1/dh3+1)
ny3 = int((ny-1)*dh1/dh3+1)
nz3 = int((nz-1)*dh1/dh3+1)

nxx1 = nx1 + 2*nb
nyy1 = ny1 + 2*nb
nzz1 = nz1 + 2*nb

nxx2 = nx2 + 2*nb
nyy2 = ny2 + 2*nb
nzz2 = nz2 + 2*nb

nxx3 = nx3 + 2*nb
nyy3 = ny3 + 2*nb
nzz3 = nz3 + 2*nb

data1 = pyf.read_binary_matrix(nrx,nry,f"eikonal_iso_{nzz1}_nStations519841_shot_1.bin")
data2 = pyf.read_binary_matrix(nrx,nry,f"eikonal_iso_{nzz2}_nStations519841_shot_1.bin")
data3 = pyf.read_binary_matrix(nrx,nry,f"eikonal_iso_{nzz3}_nStations519841_shot_1.bin")

loc = np.linspace(0, 720, 7, dtype = int)
lab = np.linspace(1, 19, 7, dtype=int)

fig, ax = plt.subplots(figsize = (10, 5), nrows = 2, ncols = 3)

im = ax[0,0].imshow(data1, cmap = "jet", vmin = np.min(data1), vmax = np.max(data1))
cbar = plt.colorbar(im, ax = ax[0,0])
cbar.set_label("Time [s]", fontsize = 15)

ax[0,0].set_title("dh = 25 m", fontsize = 18)
ax[0,0].set_xlabel("X [km]", fontsize = 15)
ax[0,0].set_ylabel("Y [km]", fontsize = 15)

ax[0,0].set_xticks(loc)
ax[0,0].set_xticklabels(lab)

ax[0,0].set_yticks(loc)
ax[0,0].set_yticklabels(lab)

ax[0,0].invert_yaxis()

im = ax[0,1].imshow(data2, cmap = "jet", vmin = np.min(data1), vmax = np.max(data1))
cbar = plt.colorbar(im, ax = ax[0,1])
cbar.set_label("Time [s]", fontsize = 15)

ax[0,1].set_title("dh = 50 m", fontsize = 18)

ax[0,1].set_xlabel("X [km]", fontsize = 15)
ax[0,1].set_ylabel("Y [km]", fontsize = 15)

ax[0,1].set_xticks(loc)
ax[0,1].set_xticklabels(lab)

ax[0,1].set_yticks(loc)
ax[0,1].set_yticklabels(lab)

ax[0,1].invert_yaxis()


im = ax[0,2].imshow(abs(data1 - data2), cmap = "jet", vmin = 0, vmax = 0.005)
cbar = plt.colorbar(im, ax = ax[0,2])
cbar.set_label("ABS(DIFF) [s]", fontsize = 15)

ax[0,2].set_title("Difference", fontsize = 18)

ax[0,2].set_xlabel("X [km]", fontsize = 15)
ax[0,2].set_ylabel("Y [km]", fontsize = 15)

ax[0,2].set_xticks(loc)
ax[0,2].set_xticklabels(lab)

ax[0,2].set_yticks(loc)
ax[0,2].set_yticklabels(lab)

ax[0,2].invert_yaxis()



im = ax[1,0].imshow(data1, cmap = "jet", vmin = np.min(data1), vmax = np.max(data1))
cbar = plt.colorbar(im, ax = ax[1,0])
cbar.set_label("Time [s]", fontsize = 15)

ax[1,0].set_title("dh = 25 m", fontsize = 18)
ax[1,0].set_xlabel("X [km]", fontsize = 15)
ax[1,0].set_ylabel("Y [km]", fontsize = 15)

ax[1,0].set_xticks(loc)
ax[1,0].set_xticklabels(lab)

ax[1,0].set_yticks(loc)
ax[1,0].set_yticklabels(lab)

ax[1,0].invert_yaxis()

im = ax[1,1].imshow(data3, cmap = "jet", vmin = np.min(data1), vmax = np.max(data1))
cbar = plt.colorbar(im, ax = ax[1,1])
cbar.set_label("Time [s]", fontsize = 15)

ax[1,1].set_title("dh = 100 m", fontsize = 18)

ax[1,1].set_xlabel("X [km]", fontsize = 15)
ax[1,1].set_ylabel("Y [km]", fontsize = 15)

ax[1,1].set_xticks(loc)
ax[1,1].set_xticklabels(lab)

ax[1,1].set_yticks(loc)
ax[1,1].set_yticklabels(lab)

ax[1,1].invert_yaxis()

im = ax[1,2].imshow(abs(data1 - data3), cmap = "jet", vmin = 0, vmax = 0.016)
cbar = plt.colorbar(im, ax = ax[1,2])
cbar.set_label("ABS(DIFF) [s]", fontsize = 15)

ax[1,2].set_title("Difference", fontsize = 18)

ax[1,2].set_xlabel("X [km]", fontsize = 15)
ax[1,2].set_ylabel("Y [km]", fontsize = 15)

ax[1,2].set_xticks(loc)
ax[1,2].set_xticklabels(lab)

ax[1,2].set_yticks(loc)
ax[1,2].set_yticklabels(lab)

ax[1,2].invert_yaxis()



fig.tight_layout()
plt.show()





# eikonal1 = pyf.read_binary_volume(nzz1,nxx1,nyy1,f"eikonal_{nzz1}.bin")
# eikonal2 = pyf.read_binary_volume(nzz2,nxx2,nyy2,f"eikonal_{nzz2}.bin")
# eikonal3 = pyf.read_binary_volume(nzz3,nxx3,nyy3,f"eikonal_{nzz3}.bin")

# eikonal1 = eikonal1[nb:nzz1-nb,nb:nxx1-nb,nb:nyy1-nb]
# eikonal2 = eikonal2[nb:nzz2-nb,nb:nxx2-nb,nb:nyy2-nb]
# eikonal3 = eikonal3[nb:nzz3-nb,nb:nxx3-nb,nb:nyy3-nb]

# slowness1 = pyf.read_binary_volume(nzz1,nxx1,nyy1,f"slowness_{nzz1}") 
# slowness2 = pyf.read_binary_volume(nzz2,nxx2,nyy2,f"slowness_{nzz2}") 
# slowness3 = pyf.read_binary_volume(nzz3,nxx3,nyy3,f"slowness_{nzz3}") 

# slowness1 = 1/slowness1[nb:nzz1-nb,nb:nxx1-nb,nb:nyy1-nb]
# slowness2 = 1/slowness2[nb:nzz2-nb,nb:nxx2-nb,nb:nyy2-nb]
# slowness3 = 1/slowness3[nb:nzz3-nb,nb:nxx3-nb,nb:nyy3-nb]

# dh = np.array([dh1, dh1, dh1])
# slices = np.array([0.25*nz1, 0.5*ny1, 0.5*nx1], dtype = int)

# sps_path = "../../FWI3D/inputs/geometry/kdm_test_SPS.txt"
# rps_path = "../../FWI3D/inputs/geometry/kdm_test_RPS.txt"

# pyf.plot_model_3D(slowness1, dh, slices, shots = sps_path, scale = 3.2, 
#                   nodes = rps_path, adjx = 0.8, dbar = 1.55, cmap = "jet",
#                   eikonal = eikonal1, eikonal_levels = np.linspace(0.25,2.0,9),
#                   cblab = "P wave velocity [km/s]")
# plt.savefig("setup_25m.png", dpi = 200)
# plt.show()

# dh = np.array([dh2, dh2, dh2])
# slices = np.array([0.25*nz2, 0.5*ny2, 0.5*nx2], dtype = int)

# pyf.plot_model_3D(slowness2, dh, slices, shots = sps_path, scale = 3.2, 
#                   nodes = rps_path, adjx = 0.8, dbar = 1.55, cmap = "jet",
#                   eikonal = eikonal2, eikonal_levels = np.linspace(0.25,2.0,9),
#                   cblab = "P wave velocity [km/s]")
# plt.savefig("setup_50m.png", dpi = 200)
# plt.show()

# dh = np.array([dh3, dh3, dh3])
# slices = np.array([0.25*nz3, 0.5*ny3, 0.5*nx3], dtype = int)

# pyf.plot_model_3D(slowness3, dh, slices, shots = sps_path, scale = 3.2, 
#                   nodes = rps_path, adjx = 0.8, dbar = 1.55, cmap = "jet",
#                   eikonal = eikonal3, eikonal_levels = np.linspace(0.25,2.0,9),
#                   cblab = "P wave velocity [km/s]")
# plt.savefig("setup_100m.png", dpi = 200)
# plt.show()
