import numpy as np

x_max = 3000
y_max = 2000
z_max = 1000

dx = 12.5
dy = 12.5
dz = 12.5

nx = int((x_max / dx) + 1)
ny = int((y_max / dy) + 1)
nz = int((z_max / dz) + 1)

A = np.array([0.5, 1, -0.5, 1, 0.5])

xc = np.array([0.25*x_max, 0.25*x_max, 0.5*x_max, 0.75*x_max, 0.75*x_max])
yc = np.array([0.25*y_max, 0.75*y_max, 0.5*y_max, 0.25*y_max, 0.75*y_max])

sigx = np.array([0.2*x_max, 0.2*x_max, 0.1*x_max, 0.2*x_max, 0.2*x_max])
sigy = np.array([0.3*y_max, 0.3*y_max, 0.1*y_max, 0.3*y_max, 0.3*y_max])

y, x = np.meshgrid(np.arange(ny)*dy, np.arange(nx)*dx)

surface = np.zeros((nx, ny))

for i in range(len(A)):
    surface += A[i]*np.exp(-0.5*( ((x - xc[i])/sigx[i])**2 + ((y - yc[i])/sigy[i])**2))

surface = z_max - 0.5*z_max/np.max(surface)*surface

vp_model = np.zeros((nz, nx, ny)) + 1500
vs_model = np.zeros((nz, nx, ny)) 
ro_model = np.zeros((nz, nx, ny)) + 1000 

vp_model[int(250/dz):] = 1700
vs_model[int(250/dz):] = 1050
ro_model[int(250/dz):] = 2260

for j in range(nx):
    for k in range(ny):
        vp_model[int(surface[j,k]/dz):, j, k] = 2000
        vs_model[int(surface[j,k]/dz):, j, k] = 2000/1.7
        ro_model[int(surface[j,k]/dz):, j, k] = 310*2000**0.25

vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_vp_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_vs_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
ro_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_ro_model_{nz}x{nx}x{ny}_{dx:.1f}m.bin")

model_input = vp_model[::2,::2,::2]
model_input.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_vp_input_{int(0.5*nz)+1}x{int(0.5*nx)+1}x{int(0.5*ny)+1}_{2.0*dx:.0f}m.bin")
