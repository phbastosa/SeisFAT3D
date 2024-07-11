import numpy as np

n = 201

homogeneous_vp_model = 1500.0*np.ones((n, n, n))
homogeneous_vs_model = 0.0*np.ones((n, n, n))
homogeneous_rho_model = 1000.0*np.ones((n, n, n))

diffraction_vp_model = homogeneous_vp_model.copy()
diffraction_vs_model = homogeneous_vs_model.copy()
diffraction_rho_model = homogeneous_rho_model.copy()

point = slice(int(0.5*n)-1,int(0.5*n)+2)

diffraction_vp_model[point,point,point] = 2000.0
diffraction_vs_model[point,point,point] = 1400.0
diffraction_rho_model[point,point,point] = 2073.0

homogeneous_vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_vp_{n}x{n}x{n}_10m.bin")
homogeneous_vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_vs_{n}x{n}x{n}_10m.bin")
homogeneous_rho_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/homogeneous_rho_{n}x{n}x{n}_10m.bin")

diffraction_vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/diffraction_vp_{n}x{n}x{n}_10m.bin")
diffraction_vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/diffraction_vs_{n}x{n}x{n}_10m.bin")
diffraction_rho_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/diffraction_rho_{n}x{n}x{n}_10m.bin")
