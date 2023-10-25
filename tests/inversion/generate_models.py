import numpy as np

x = 5000
y = 5000
z = 2000

dh = 50.0

nx = int((x / dh) + 1)
ny = int((y / dh) + 1)
nz = int((z / dh) + 1)

init_model = 1500.0 * np.ones((nz,nx,ny))
true_model = 1500.0 * np.ones((nz,nx,ny))

zslice = slice(int(750/dh),int(1250/dh))
xslice = slice(int(1500/dh),int(3500/dh))
yslice = slice(int(1500/dh),int(3500/dh))

true_model[zslice,xslice,yslice] = 2000.0

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")