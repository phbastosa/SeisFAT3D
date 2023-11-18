import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------------

def createGaussianSurface(nx,ny,dh,A,xc,yc,sigx,sigy):

    x, y = np.meshgrid(np.arange(nx)*dh, np.arange(ny)*dh)

    surface = A*np.exp(-((x - xc)**2/(2*sigx**2) + (y - yc)**2/(2*sigy**2)))

    return surface

#------------------------------------------------------------------------------------

dh = 10

nx = 701
ny = 501
nz = 151

nlayers = 11
initial = 1500.0
increase = 100.0

basinDepth = 1300.0
waterBottom = 300.0

surface_top = 400.0
surface_base = 1400.0
surface_velocity = 3000.0

A = np.array([650, 800, -550, 900, 650])

xc = np.array([2500, 2500, 3500, 4500, 4500])
yc = np.array([2000, 3000, 2500, 2000, 3000])

sigx = np.array([800, 800, 500, 800, 800])
sigy = np.array([500, 500, 300, 600, 600])

#------------------------------------------------------------------------------------

velocities = initial + np.arange(nlayers)*increase 

interfaces = np.linspace(waterBottom, basinDepth, nlayers-1)

surface = np.zeros((ny, nx))

model = initial * np.ones((nz, nx, ny))

for index, depth in enumerate(interfaces):
    model[int(depth/dh):] = velocities[index+1]

for k in range(len(A)):
    surface += createGaussianSurface(nx, ny, dh, A[k], xc[k], yc[k], sigx[k], sigy[k])

surface = surface_top//dh + np.array((np.max(surface) - 0.8*surface)//dh, dtype=int)

i,j = np.where(-surface*dh < -surface_base)

surface[i,j] = surface_base/dh

#------------------------------------------------------------------------------------

_, depth, _ = np.meshgrid(np.arange(nx), np.arange(nz), np.arange(ny))

update = depth >= surface.T

model[update] = surface_velocity

#------------------------------------------------------------------------------------

model.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/refractive_model_{nz}x{nx}x{ny}_{dh:.0f}m.bin")