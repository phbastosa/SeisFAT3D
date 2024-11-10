import numpy as np
import matplotlib.pyplot as plt

x_max = 3e3
y_max = 2e3
z_max = 5e2

dx = 12.5
dy = 12.5
dz = 12.5

nx = int((x_max / dx) + 1)
ny = int((y_max / dy) + 1)
nz = int((z_max / dz) + 1)

A = np.array([1, 2, -1, 2, 1])

xc = np.array([0.25*x_max, 0.25*x_max, 0.5*x_max, 0.75*x_max, 0.75*x_max])
yc = np.array([0.25*y_max, 0.75*y_max, 0.5*y_max, 0.25*y_max, 0.75*y_max])

sigx = np.array([0.2*x_max, 0.2*x_max, 0.1*x_max, 0.2*x_max, 0.2*x_max])
sigy = np.array([0.3*y_max, 0.3*y_max, 0.1*y_max, 0.3*y_max, 0.3*y_max])

y, x = np.meshgrid(np.arange(ny)*dy, np.arange(nx)*dx)

surface = np.zeros((nx, ny))

for i in range(len(A)):
    surface += A[i]*np.exp(-0.5*( ((x - xc[i])/sigx[i])**2 + ((y - yc[i])/sigy[i])**2))

surface = z_max - 0.75*z_max/np.max(surface)*surface

vp_model = np.zeros((nz, nx, ny)) + 2000

for j in range(nx):
    for k in range(ny):
        vp_model[int(surface[j,k]/dz):, j, k] = 2500.0

vs_model = 0.7*vp_model
rho_model = 310*vp_model**0.25

vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_vp_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_vs_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
rho_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_rho_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
