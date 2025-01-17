import numpy as np

nx = 201
ny = 101
nz = 51

dh = 10.0

# Thomsen paramters for Orthorhombic case

E1 = 0.0
E2 = 0.0
D1 = 0.0
D2 = 0.0
D3 = 0.0
G1 = 0.0
G2 = 0.0

# Tilt application

theta_x = 0.0 * np.pi/180.0
theta_y = 0.0 * np.pi/180.0
theta_z = 0.0 * np.pi/180.0

# Elastic properties

Vp = 2000.0
Vs = 1100.0
Ro = 2600.0

C11 = C12 = C13 = C14 = C15 = C16 = 0
C22 = C23 = C24 = C25 = C26 = 0
C33 = C34 = C35 = C36 = 0
C44 = C45 = C46 = 0
C55 = C56 = 0
C66 = 0

# Common constants 

C33 = Ro*Vp**2
C55 = Ro*Vs**2

# Orthorhombic case

C11 = C33*(1.0 + 2.0*E2)
C22 = C33*(1.0 + 2.0*E1)

C66 = C55*(1.0 + 2.0*G1)
C44 = C66/(1.0 + 2.0*G2)

C23 = np.sqrt((C33 - C44)**2 + 2.0*D1*C33*(C33 - C44)) - C44
C13 = np.sqrt((C33 - C55)**2 + 2.0*D2*C33*(C33 - C55)) - C55
C12 = np.sqrt((C11 - C66)**2 + 2.0*D3*C11*(C11 - C66)) - C66

# Tilt application

# C = np.zeros((6,6))
# M = np.zeros((6,6))

# C[0,0] = C11; C[0,1] = C12; C[0,2] = C13; C[0,3] = C14; C[0,4] = C15; C[0,5] = C16  
# C[1,0] = C12; C[1,1] = C22; C[1,2] = C23; C[1,3] = C24; C[1,4] = C25; C[1,5] = C26  
# C[2,0] = C13; C[2,1] = C12; C[2,2] = C33; C[2,3] = C34; C[2,4] = C35; C[2,5] = C36  
# C[3,0] = C14; C[3,1] = C24; C[3,2] = C34; C[3,3] = C44; C[3,4] = C45; C[3,5] = C46  
# C[4,0] = C15; C[4,1] = C25; C[4,2] = C35; C[4,3] = C45; C[4,4] = C55; C[4,5] = C56  
# C[5,0] = C16; C[5,1] = C26; C[5,2] = C36; C[5,3] = C46; C[5,4] = C56; C[5,5] = C66  

# Rx = np.array([
#     [1, 0, 0],
#     [0, np.cos(theta_x),-np.sin(theta_x)],
#     [0, np.sin(theta_x), np.cos(theta_x)]
# ])
    
# Ry = np.array([
#     [np.cos(theta_y), 0, np.sin(theta_y)],
#     [0, 1, 0],
#     [-np.sin(theta_y), 0, np.cos(theta_y)]
# ])

# Rz = np.array([
#     [np.cos(theta_z),-np.sin(theta_z), 0],
#     [np.sin(theta_z), np.cos(theta_z), 0],
#     [0, 0, 1]
# ])

# R = Rz @ Ry @ Rx

# pax = R[0,0]; pay = R[0,1]; paz = R[0,2]
# pbx = R[1,0]; pby = R[1,1]; pbz = R[1,2]
# pcx = R[2,0]; pcy = R[2,1]; pcz = R[2,2]

# M[0,0] = pax*pax; M[1,0] = pbx*pbx; M[2,0] = pcx*pcx
# M[0,1] = pay*pay; M[1,1] = pby*pby; M[2,1] = pcy*pcy
# M[0,2] = paz*paz; M[1,2] = pbz*pbz; M[2,2] = pcz*pcz

# M[3,0] = pbx*pcx; M[4,0] = pcx*pax; M[5,0] = pax*pbx
# M[3,1] = pby*pcy; M[4,1] = pcy*pay; M[5,1] = pay*pby
# M[3,2] = pbz*pcz; M[4,2] = pcz*paz; M[5,2] = paz*pbz

# M[0,3] = 2.0*pay*paz; M[1,3] = 2.0*pby*pbz; M[2,3] = 2.0*pcy*pcz
# M[0,4] = 2.0*paz*pax; M[1,4] = 2.0*pbz*pbx; M[2,4] = 2.0*pcz*pcx
# M[0,5] = 2.0*pax*pay; M[1,5] = 2.0*pbx*pby; M[2,5] = 2.0*pcx*pcy

# M[3,3] = pby*pcz + pbz*pcy; M[4,3] = pay*pcz + paz*pcy; M[5,3] = pay*pbz + paz*pby
# M[3,4] = pbx*pcz + pbz*pcx; M[4,4] = paz*pcx + pax*pcz; M[5,4] = paz*pbx + pax*pbz
# M[3,5] = pby*pcx + pbx*pcy; M[4,5] = pax*pcy + pay*pcx; M[5,5] = pax*pby + pay*pbx

# C_rot = M @ C @ M.T

C_rot = np.array([[C11,C12,C13,C14,C15,C16],
                  [ 0 ,C22,C23,C24,C25,C26],
                  [ 0 , 0 ,C33,C34,C35,C36],
                  [ 0 , 0 , 0 ,C44,C45,C46],
                  [ 0 , 0 , 0 , 0 ,C55,C56],
                  [ 0 , 0 , 0 , 0 , 0 ,C66]]) 

# Filling properties

Vp = np.zeros((nz,nx,ny)) + Vp
Vs = np.zeros((nz,nx,ny)) + Vs
Ro = np.zeros((nz,nx,ny)) + Ro

C11 = np.zeros_like(Vp) + C_rot[0,0]
C12 = np.zeros_like(Vp) + C_rot[0,1]
C13 = np.zeros_like(Vp) + C_rot[0,2]
C14 = np.zeros_like(Vp) + C_rot[0,3]
C15 = np.zeros_like(Vp) + C_rot[0,4]
C16 = np.zeros_like(Vp) + C_rot[0,5]

C22 = np.zeros_like(Vp) + C_rot[1,1]
C23 = np.zeros_like(Vp) + C_rot[1,2]
C24 = np.zeros_like(Vp) + C_rot[1,3]
C25 = np.zeros_like(Vp) + C_rot[1,4]
C26 = np.zeros_like(Vp) + C_rot[1,5]

C33 = np.zeros_like(Vp) + C_rot[2,2]
C34 = np.zeros_like(Vp) + C_rot[2,3]
C35 = np.zeros_like(Vp) + C_rot[2,4]
C36 = np.zeros_like(Vp) + C_rot[2,5]

C44 = np.zeros_like(Vp) + C_rot[3,3]
C45 = np.zeros_like(Vp) + C_rot[3,4]
C46 = np.zeros_like(Vp) + C_rot[3,5]

C55 = np.zeros_like(Vp) + C_rot[4,4]
C56 = np.zeros_like(Vp) + C_rot[4,5]

C66 = np.zeros_like(Vp) + C_rot[5,5]

# Export models

Vp.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_vp.bin")
Vs.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_vs.bin")
Ro.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_ro.bin")

C11.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C11.bin")
C12.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C12.bin")
C13.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C13.bin")
C14.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C14.bin")
C15.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C15.bin")
C16.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C16.bin")

C22.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C22.bin")
C23.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C23.bin")
C24.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C24.bin")
C25.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C25.bin")
C26.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C26.bin")

C33.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C33.bin")
C34.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C34.bin")
C35.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C35.bin")
C36.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C36.bin")

C44.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C44.bin")
C45.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C45.bin")
C46.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C46.bin")

C55.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C55.bin")
C56.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C56.bin")

C66.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTest_C66.bin")

# Acquisition geometry

nsx = 1
nsy = 1

nrx = 181 
nry = 5

ns = nsx*nsy
nr = nrx*nry

sx, sy = 50, 500 
rx, ry = np.meshgrid(np.linspace(100, 1900, nrx), np.linspace(100, 900, nry))

SPS = np.zeros((ns, 3), dtype = float)
RPS = np.zeros((nr, 3), dtype = float)
XPS = np.zeros((ns, 3), dtype = int)

SPS[:,0] = np.reshape(sx, [ns], order = "F")
SPS[:,1] = np.reshape(sy, [ns], order = "F")
SPS[:,2] = np.zeros(ns) + dh

RPS[:,0] = np.reshape(rx, [nr], order = "C")
RPS[:,1] = np.reshape(ry, [nr], order = "C")
RPS[:,2] = np.zeros(nr) + dh*(nz-1) - dh

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.zeros(ns) 
XPS[:, 2] = np.zeros(ns) + nr 

np.savetxt("../inputs/geometry/anisoTest_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTest_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTest_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")
