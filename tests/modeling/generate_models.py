import numpy as np

x_max = 8e3
y_max = 5e3
z_max = 2e3

dh = 25.0

nx = int((x_max / dh) + 1)
ny = int((y_max / dh) + 1)
nz = int((z_max / dh) + 1)

vp_model = np.zeros((nz, nx, ny)) + 1500
vs_model = np.zeros((nz, nx, ny)) 
rho_model = np.zeros((nz, nx, ny)) + 1000

v = np.array([1500, 1700, 1900, 2300, 3000])
z = np.array([400, 400, 400, 400])

for i in range(len(z)):
    vp_model[int(np.sum(z[:i+1]/dh)):] = v[i+1]
    vs_model[int(np.sum(z[:i+1]/dh)):] = 0.7*v[i+1]
    rho_model[int(np.sum(z[:i+1]/dh)):] = 310*v[i+1]**0.25

vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vp_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vs_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
rho_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_rho_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
