import numpy as np
import matplotlib.pyplot as plt

from sys import path, argv
path.append("../src/")
import functions

nx = int(functions.catch_parameter(argv[1], "x_samples"))
ny = int(functions.catch_parameter(argv[1], "y_samples"))
nz = int(functions.catch_parameter(argv[1], "z_samples"))

dx = float(functions.catch_parameter(argv[1], "x_spacing"))
dy = float(functions.catch_parameter(argv[1], "y_spacing"))
dz = float(functions.catch_parameter(argv[1], "z_spacing"))

shots_file = functions.catch_parameter(argv[1], "shots_file")
nodes_file = functions.catch_parameter(argv[1], "nodes_file")

vp_model_file = functions.catch_parameter(argv[1], "vp_model_file")
vs_model_file = functions.catch_parameter(argv[1], "vs_model_file")
rho_model_file = functions.catch_parameter(argv[1], "rho_model_file")

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
                        scale = 2.0)

plt.savefig(f"vp_model_test.png", dpi = 200)
plt.clf()

functions.plot_model_3D(model_vs, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 3000,
                        scale = 2.0)

plt.savefig(f"vs_model_test.png", dpi = 200)
plt.clf()

functions.plot_model_3D(model_rho, dh, slices,
                        shots = shots_file,
                        nodes = nodes_file,
                        vmin = 0,
                        vmax = 3000,
                        scale = 2.0)

plt.savefig(f"rho_model_test.png", dpi = 200)
plt.clf()

nTraces = len(nodes)

nt = int(functions.catch_parameter(argv[1], "time_samples"))
dt = float(functions.catch_parameter(argv[1], "time_spacing"))

scalar_seismogram_path = f"../outputs/seismograms/scalar_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"
elastic_seismogram_path = f"../outputs/seismograms/elastic_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"
acoustic_seismogram_path = f"../outputs/seismograms/acoustic_seismogram_Nsamples{nt}_nRec{nTraces}_shot_1.bin"

scalar_seism = functions.read_binary_matrix(nt, nTraces, scalar_seismogram_path)
elastic_seism = functions.read_binary_matrix(nt, nTraces, elastic_seismogram_path)
acoustic_seism = functions.read_binary_matrix(nt, nTraces, acoustic_seismogram_path)

scale = 10.0*np.std(elastic_seism)

fig, ax = plt.subplots(figsize = (10,6))

plt.imshow(scalar_seism, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)

fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize = (10,6))

plt.imshow(acoustic_seism, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)

fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize = (10,6))

plt.imshow(elastic_seism, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)

fig.tight_layout()
plt.show()
