import numpy as np

# Geometry configurations ------------------------------------------------------------------

xc = 10000
yc = 10000
ds = 100.0

offsets = [7000, 8000, 9000]

sx = np.array([])
sy = np.array([])

for offset in offsets:

    angles = np.linspace(0.0, 2.0*np.pi, int(2.0*np.pi*offset/ds))

    sx = np.append(sx, offset * np.cos(angles) + xc)
    sy = np.append(sy, offset * np.sin(angles) + yc)

nrx = 19
nry = 19

rx, ry = np.meshgrid(np.linspace(5500, 14500, nrx), np.linspace(5500, 14500, nrx))

# Writing with reciprocity principle

ns = nrx*nry
nr = len(sx)

SPS = np.zeros((ns, 3), dtype = float)
RPS = np.zeros((nr, 3), dtype = float)
XPS = np.zeros((ns, 3), dtype = int)

SPS[:,0] = np.reshape(rx, [ns], order = "F")
SPS[:,1] = np.reshape(ry, [ns], order = "F")
SPS[:,2] = 1000.0

RPS[:,0] = sx
RPS[:,1] = sy
RPS[:,2] = 8.0

np.savetxt("../inputs/geometry/anisoTomo_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTomo_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")

# Model configurations ---------------------------------------------------------------------

x_max = 2e4
y_max = 2e4
z_max = 5e3

dh = 100.0

nx = int((x_max / dh) + 1)
ny = int((y_max / dh) + 1)
nz = int((z_max / dh) + 1)

vp = np.zeros((nz, nx, ny)) + 1500
vs = np.zeros((nz, nx, ny))
ro = np.zeros((nz, nx, ny)) + 1000

E1 = np.zeros((nz, nx, ny))
E2 = np.zeros((nz, nx, ny))

D1 = np.zeros((nz, nx, ny))
D2 = np.zeros((nz, nx, ny))
D3 = np.zeros((nz, nx, ny))

G1 = np.zeros((nz, nx, ny))
G2 = np.zeros((nz, nx, ny))

dip = np.zeros((nz, nx, ny))
azm = np.zeros((nz, nx, ny))

dv = 50.0
vi = 1650.0
wb = 1000.0

for i in range(nz):
    if i > wb/dh:
        vp[i] = vi + (i*dh - wb)*dv/dh 
        vs[i] = vp[i] / 1.7
        ro[i] = 310*ro[i]**0.25    

vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/anisoTomo_vp_init.bin")

radius = 1000

dvp = np.array([500, -500, -500, 500])
dvs = np.array([300, -300, -300, 300])
dro = np.array([ 50,  -50,  -50,  50])

dE1 = np.array([0.0, 0.0, 0.0, 0.0])
dE2 = np.array([0.0, 0.0, 0.0, 0.0])

dD1 = np.array([0.0, 0.0, 0.0, 0.0])
dD2 = np.array([0.0, 0.0, 0.0, 0.0])
dD3 = np.array([0.0, 0.0, 0.0, 0.0])

dG1 = np.array([0.0, 0.0, 0.0, 0.0])
dG2 = np.array([0.0, 0.0, 0.0, 0.0])

ddip = np.deg2rad(np.array([0.0, 0.0, 0.0, 0.0]))
dazm = np.deg2rad(np.array([0.0, 0.0, 0.0, 0.0]))

circle_centers = np.array([[3000, 8000, 8000],
                           [3000, 8000, 12000],
                           [3000, 12000, 8000],
                           [3000, 12000, 12000]])

x, z, y = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh, np.arange(ny)*dh)

points = []

for k in range(len(circle_centers)):
    
    distance = np.sqrt((z - circle_centers[k,0])**2 + (x - circle_centers[k,1])**2 + (y - circle_centers[k,2])**2)

    vp[distance <= radius] += dvp[k]
    vs[distance <= radius] += dvs[k]
    vs[distance <= radius] += dvs[k]

    E1[distance <= radius] += dE1[k]
    E2[distance <= radius] += dE2[k]

    D1[distance <= radius] += dD1[k]
    D2[distance <= radius] += dD2[k]    
    D3[distance <= radius] += dD3[k]

    G1[distance <= radius] += dG1[k]    
    G2[distance <= radius] += dG2[k]

    dip[distance <= radius] += ddip[k]
    azm[distance <= radius] += dazm[k]

    points.append(np.where(distance <= radius))

C11 = np.zeros_like(vp)
C12 = np.zeros_like(vp) 
C13 = np.zeros_like(vp) 
C14 = np.zeros_like(vp) 
C15 = np.zeros_like(vp)
C16 = np.zeros_like(vp)

C22 = np.zeros_like(vp)
C23 = np.zeros_like(vp)
C24 = np.zeros_like(vp)
C25 = np.zeros_like(vp)
C26 = np.zeros_like(vp)

C33 = np.zeros_like(vp)
C34 = np.zeros_like(vp)
C35 = np.zeros_like(vp)
C36 = np.zeros_like(vp)

C44 = np.zeros_like(vp)
C45 = np.zeros_like(vp)
C46 = np.zeros_like(vp)

C55 = np.zeros_like(vp)
C56 = np.zeros_like(vp)

C66 = np.zeros_like(vp)

C = np.zeros((6,6))
M = np.zeros((6,6))

c11 = c12 = c13 = c14 = c15 = c16 = 0
c22 = c23 = c24 = c25 = c26 = 0
c33 = c34 = c35 = c36 = 0
c44 = c45 = c46 = 0
c55 = c56 = 0
c66 = 0

SI = 1e9

C33 = ro*vp**2 / SI
C55 = ro*vs**2 / SI

C11 = C33*(1.0 + 2.0*E2)
C22 = C33*(1.0 + 2.0*E1)

C66 = C55*(1.0 + 2.0*G1)
C44 = C66/(1.0 + 2.0*G2)

C23 = np.sqrt((C33 - C44)**2 + 2.0*D1*C33*(C33 - C44)) - C44
C13 = np.sqrt((C33 - C55)**2 + 2.0*D2*C33*(C33 - C55)) - C55
C12 = np.sqrt((C11 - C66)**2 + 2.0*D3*C11*(C11 - C66)) - C66

for circles in points:
    
    ip, jp, kp = circles
    
    for index in range(len(ip)):
        
        i = ip[index]
        j = jp[index]
        k = kp[index]

        C[0,0] = C11[i,j,k]; C[0,1] = C12[i,j,k]; C[0,2] = C13[i,j,k]; C[0,3] = C14[i,j,k]; C[0,4] = C15[i,j,k]; C[0,5] = C16[i,j,k]  
        C[1,0] = C12[i,j,k]; C[1,1] = C22[i,j,k]; C[1,2] = C23[i,j,k]; C[1,3] = C24[i,j,k]; C[1,4] = C25[i,j,k]; C[1,5] = C26[i,j,k]  
        C[2,0] = C13[i,j,k]; C[2,1] = C23[i,j,k]; C[2,2] = C33[i,j,k]; C[2,3] = C34[i,j,k]; C[2,4] = C35[i,j,k]; C[2,5] = C36[i,j,k]  
        C[3,0] = C14[i,j,k]; C[3,1] = C24[i,j,k]; C[3,2] = C34[i,j,k]; C[3,3] = C44[i,j,k]; C[3,4] = C45[i,j,k]; C[3,5] = C46[i,j,k]  
        C[4,0] = C15[i,j,k]; C[4,1] = C25[i,j,k]; C[4,2] = C35[i,j,k]; C[4,3] = C45[i,j,k]; C[4,4] = C55[i,j,k]; C[4,5] = C56[i,j,k]  
        C[5,0] = C16[i,j,k]; C[5,1] = C26[i,j,k]; C[5,2] = C36[i,j,k]; C[5,3] = C46[i,j,k]; C[5,4] = C56[i,j,k]; C[5,5] = C66[i,j,k]  

        c = np.cos(dip[i,j,k]); s = np.sin(dip[i,j,k])
        Rdip = np.array([[c, 0, s],[0, 1, 0],[-s, 0, c]])

        c = np.cos(azm[i,j,k]); s = np.sin(azm[i,j,k])
        Razm = np.array([[c,-s, 0],[s, c, 0],[0, 0, 1]])

        R = Rdip @ Razm

        pax = R[0,0]; pay = R[0,1]; paz = R[0,2]
        pbx = R[1,0]; pby = R[1,1]; pbz = R[1,2]
        pcx = R[2,0]; pcy = R[2,1]; pcz = R[2,2]

        M[0,0] = pax*pax; M[1,0] = pbx*pbx; M[2,0] = pcx*pcx
        M[0,1] = pay*pay; M[1,1] = pby*pby; M[2,1] = pcy*pcy
        M[0,2] = paz*paz; M[1,2] = pbz*pbz; M[2,2] = pcz*pcz

        M[3,0] = pbx*pcx; M[4,0] = pcx*pax; M[5,0] = pax*pbx
        M[3,1] = pby*pcy; M[4,1] = pcy*pay; M[5,1] = pay*pby
        M[3,2] = pbz*pcz; M[4,2] = pcz*paz; M[5,2] = paz*pbz

        M[0,3] = 2.0*pay*paz; M[1,3] = 2.0*pby*pbz; M[2,3] = 2.0*pcy*pcz
        M[0,4] = 2.0*paz*pax; M[1,4] = 2.0*pbz*pbx; M[2,4] = 2.0*pcz*pcx
        M[0,5] = 2.0*pax*pay; M[1,5] = 2.0*pbx*pby; M[2,5] = 2.0*pcx*pcy

        M[3,3] = pby*pcz + pbz*pcy; M[4,3] = pay*pcz + paz*pcy; M[5,3] = pay*pbz + paz*pby
        M[3,4] = pbx*pcz + pbz*pcx; M[4,4] = paz*pcx + pax*pcz; M[5,4] = paz*pbx + pax*pbz
        M[3,5] = pby*pcx + pbx*pcy; M[4,5] = pax*pcy + pay*pcx; M[5,5] = pax*pby + pay*pbx

        Cr = (M @ C @ M.T) * SI

        C11[i,j,k] = Cr[0,0]; C12[i,j,k] = Cr[0,1]; C13[i,j,k] = Cr[0,2]; C14[i,j,k] = Cr[0,3]
        C15[i,j,k] = Cr[0,4]; C16[i,j,k] = Cr[0,5]; C22[i,j,k] = Cr[1,1]; C23[i,j,k] = Cr[1,2]
        C24[i,j,k] = Cr[1,3]; C25[i,j,k] = Cr[1,4]; C26[i,j,k] = Cr[1,5]; C33[i,j,k] = Cr[2,2]
        C34[i,j,k] = Cr[2,3]; C35[i,j,k] = Cr[2,4]; C36[i,j,k] = Cr[2,5]; C44[i,j,k] = Cr[3,3]
        C45[i,j,k] = Cr[3,4]; C46[i,j,k] = Cr[3,5]; C55[i,j,k] = Cr[4,4]; C56[i,j,k] = Cr[4,5]
        C66[i,j,k] = Cr[5,5]; 

vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/anisoTomo_vp_true.bin")

C11.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C11.bin")
C12.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C12.bin")
C13.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C13.bin")
C14.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C14.bin")
C15.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C15.bin")
C16.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C16.bin")

C22.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C22.bin")
C23.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C23.bin")
C24.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C24.bin")
C25.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C25.bin")
C26.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C26.bin")

C33.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C33.bin")
C34.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C34.bin")
C35.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C35.bin")
C36.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C36.bin")

C44.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C44.bin")
C45.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C45.bin")
C46.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C46.bin")

C55.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C55.bin")
C56.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C56.bin")

C66.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/anisoTomo_C66.bin")