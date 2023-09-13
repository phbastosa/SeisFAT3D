import numpy as np

nx = 201
ny = 201
nz = 51

dh = 10

true_model = 1500.0 * np.ones((nz, nx, ny))
init_model = 1500.0 * np.ones((nz, nx, ny))

true_model[20:31,50:151,50:151] = 2000.0

true_model.flatten("F").astype("float32", order = "F").tofile(f"true_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
init_model.flatten("F").astype("float32", order = "F").tofile(f"init_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")