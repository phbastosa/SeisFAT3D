import numpy as np

nx = 301
ny = 201
nz = 101

dh = 10.0

v = np.array([1500, 1750, 2000], dtype = float)
z = np.array((400, 800), dtype = float)

model_vp = np.ones((nz,nx,ny)) * 1500.0 
model_vs = np.ones((nz,nx,ny)) * 0.0000
model_pb = np.ones((nz,nx,ny)) * 1000.0

for i in range(len(z)):
    model_vp[int(z[i]/dh):] = v[i+1]
    model_vs[int(z[i]/dh):] = v[i+1] / 1.7
    model_pb[int(z[i]/dh):] = 310.0*v[i+1]**0.23  
    
model_vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/layercake_vp_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
model_vs.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/layercake_vs_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
model_pb.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/layercake_rho_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

