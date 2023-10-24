import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def readBinaryArray(n,filename):
    return np.fromfile(filename, dtype = np.float32, count = n)

def readBinaryMatrix(n1,n2,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2)    
    return np.reshape(data, [n1,n2], order='F')

def readBinaryVolume(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

def analytical_firstArrival(v, z, x):
    direct_wave = x / v[0]

    first_arrivals = np.zeros(len(x))
    refracted_waves = np.zeros((len(z), len(x)))

    for n in range(len(z)):
        refracted_waves[n,:] = x / v[n+1]
        for i in range(n+1):
            angle = np.arcsin(v[i] / v[n+1])
            refracted_waves[n,:] += 2.0*z[i]*np.cos(angle) / v[i]
    
    for offset in range(len(x)):
        first_arrivals[offset] = np.min(np.append(refracted_waves[:, offset], direct_wave[offset]))

    return first_arrivals

def check_geometry(models, shots, nodes, dh, slices, subplots, scale = 2.0):
    
    if np.sum(subplots) == 2:
        modelShape = np.array(np.shape(models))
        maxModelDistance = np.max(np.shape(models))
        minModelDistance = np.min(np.shape(models))

        vmin = np.min(models)
        vmax = np.max(models)

    else:
        modelShape = np.array(np.shape(models[0]))
        maxModelDistance = np.max(np.shape(models[0]))
        minModelDistance = np.min(np.shape(models[0]))

        vmin = np.min(models[0])
        vmax = np.max(models[0])

    nz, nx, ny = modelShape
    [z, x, y] = scale * (minModelDistance / maxModelDistance) * modelShape / maxModelDistance

    px = 1/plt.rcParams['figure.dpi']  
    ticks = np.array([3,7,7], dtype = int)

    fig = plt.figure(1, figsize=(910*px*subplots[1], 780*px*subplots[0]))

    xloc = np.linspace(0,nx-1,ticks[1], dtype = int)
    yloc = np.linspace(0,ny-1,ticks[2], dtype = int)
    zloc = np.linspace(0,nz-1,ticks[0], dtype = int)

    m2km = 1e-3

    xlab = np.around(xloc * dh[0] * m2km, decimals = 1)
    ylab = np.around(yloc * dh[1] * m2km, decimals = 1)
    zlab = np.around(zloc * dh[2] * m2km, decimals = 1)

    axes = np.array([[0.75 - x, 0.98 - y      , x, y], 
                     [    0.75, 0.98 - y      , z, y],
                     [0.75 - x, 0.98 - y - z  , x, z],
                     [0.75 - x, 0.98 - y - 1.8*z, x, z]])

    xTickDirection = ['out', 'out', 'out']
    yTickDirection = ['out', 'in', 'out']

    xTickLock = [xloc, zloc[1:], xloc]
    yTickLock = [yloc, yloc, zloc[1:]]

    xTickLabel = [[], zlab[1:], xlab]
    yTickLabel = [ylab, [], zlab[1:]]

    xLabel = ["X [km]", "Z [km]", "X [km]"]
    yLabel = ["Y [km]", "      ", "Z [km]"]

    yInvert = [True, True, False]

    xSlices = [[np.arange(modelShape[1]), np.ones(modelShape[1])*slices[1], "--g"],
               [np.arange(modelShape[0]), np.ones(modelShape[0])*slices[1], "--g"],
               [np.arange(modelShape[1]), np.ones(modelShape[1])*slices[0], "--r"]] 

    ySlices = [[np.ones(modelShape[2])*slices[2], np.arange(modelShape[2]), "--m"],
               [np.ones(modelShape[2])*slices[0], np.arange(modelShape[2]), "--r"],
               [np.ones(modelShape[0])*slices[2], np.arange(modelShape[0]), "--m"]]

    # picking geometry     

    zy_plane_shot_y = np.array([])
    zy_plane_shot_z = np.array([])

    if np.size(shots) != 3:
        for i in range(len(shots)):    
            if int(shots[i,0]/dh[0])-1 <= int(slices[2]) <= int(shots[i,0]/dh[0])+1:
                zy_plane_shot_y = np.append(zy_plane_shot_y, shots[i,1]/dh[1])        
                zy_plane_shot_z = np.append(zy_plane_shot_z, shots[i,2]/dh[2])        
    else:
        if int(shots[0]/dh[0])-1 <= int(slices[2]) <= int(shots[0]/dh[0])+1:
            zy_plane_shot_y = np.append(zy_plane_shot_y, shots[1]/dh[1])        
            zy_plane_shot_z = np.append(zy_plane_shot_z, shots[2]/dh[2])        

    zx_plane_shot_x = np.array([])
    zx_plane_shot_z = np.array([]) 

    if np.size(shots) != 3:
        for i in range(len(shots)):
            if int(shots[i,1]/dh[1])-1 <= int(slices[1]) <= int(shots[i,1]/dh[1])+1:
                zx_plane_shot_x = np.append(zx_plane_shot_x, shots[i,0]/dh[0])        
                zx_plane_shot_z = np.append(zx_plane_shot_z, shots[i,2]/dh[2])        
    else:
        if int(shots[1]/dh[1])-1 <= int(slices[1]) <= int(shots[1]/dh[1])+1:
            zx_plane_shot_x = np.append(zx_plane_shot_x, shots[0]/dh[0])        
            zx_plane_shot_z = np.append(zx_plane_shot_z, shots[2]/dh[2])        
        
    zy_plane_node_y = np.array([])
    zy_plane_node_z = np.array([])

    if np.size(nodes) != 3:
        for i in range(len(nodes)):
            if int(nodes[i,0]/dh[0])-1 <= int(slices[2]) <= int(nodes[i,0]/dh[0])+1:
                zy_plane_node_y = np.append(zy_plane_node_y, nodes[i,1]/dh[1])        
                zy_plane_node_z = np.append(zy_plane_node_z, nodes[i,2]/dh[2])        
    else:
        if int(nodes[0]/dh[0])-1 <= int(slices[2]) <= int(nodes[0]/dh[0])+1:
            zy_plane_node_y = np.append(zy_plane_node_y, nodes[1]/dh[1])        
            zy_plane_node_z = np.append(zy_plane_node_z, nodes[2]/dh[2])        

    zx_plane_node_x = np.array([])
    zx_plane_node_z = np.array([]) 

    if np.size(nodes) != 3:
        for i in range(len(nodes)):
            if int(nodes[i,1]/dh[1])-1 <= int(slices[1]) <= int(nodes[i,1]/dh[1])+1:
                zx_plane_node_x = np.append(zx_plane_node_x, nodes[i,0]/dh[0])        
                zx_plane_node_z = np.append(zx_plane_node_z, nodes[i,2]/dh[2])        
    else:
        if int(nodes[1]/dh[1])-1 <= int(slices[1]) <= int(nodes[1]/dh[1])+1:
            zx_plane_node_x = np.append(zx_plane_node_x, nodes[0]/dh[0])        
            zx_plane_node_z = np.append(zx_plane_node_z, nodes[2]/dh[2])        
        
    #--------------------------------------------------------------------------------    

    subfigs = fig.subfigures(subplots[0], subplots[1])
    
    for i in range(subplots[0]):
        for j in range(subplots[1]):

            ind = i*subplots[0] + j 

            if np.sum(subplots) == 2:
                ims = [models[slices[0],:,:].T, models[:,slices[2],:].T, models[:,:,slices[1]]]
            else:
                ims = [models[ind, slices[0],:,:].T, models[ind,:,slices[2],:].T, models[ind,:,:,slices[1]]]

            if np.size(shots) != 3:
                xshot = [shots[:,0]/dh[0],zy_plane_shot_z,zx_plane_shot_x]
                yshot = [shots[:,1]/dh[1],zy_plane_shot_y,zx_plane_shot_z]
            else:
                xshot = [shots[0]/dh[0],zy_plane_shot_z,zx_plane_shot_x]
                yshot = [shots[1]/dh[1],zy_plane_shot_y,zx_plane_shot_z]

            if np.size(nodes) != 3:
                xnode = [nodes[:,0]/dh[0],zy_plane_node_z,zx_plane_node_x]
                ynode = [nodes[:,1]/dh[1],zy_plane_node_y,zx_plane_node_z]
            else:
                xnode = [nodes[0]/dh[0],zy_plane_node_z,zx_plane_node_x]
                ynode = [nodes[1]/dh[1],zy_plane_node_y,zx_plane_node_z]

            for k, axs in enumerate(axes):      
                if subplots[0] == 1:
                    if subplots[1] == 1:
                        ax = subfigs.add_axes(axs)                         
                    else:
                        ax = subfigs[j].add_axes(axs)
                elif subplots[1] == 1:
                    if subplots[0] == 1:
                        ax = subfigs.add_axes(axs)        
                    else:    
                        ax = subfigs[i].add_axes(axs)
                else:
                    ax = subfigs[i,j].add_axes(axs)

                # Setting colorbar
                if k == 3:
                    ax.axis("off")
                    cmap = mpl.colormaps["jet"]
                    norm = mpl.colors.Normalize(vmin*1e-3, vmax*1e-3)
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("bottom", size="10%", pad=0)
                    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(vmin*1e-3, vmax*1e-3, 5), orientation = "horizontal")
                    cbar.ax.set_xticklabels(np.around(np.linspace(vmin*1e-3, vmax*1e-3, 5), decimals = 1))
                    cbar.set_label("Velocity [km/s]")
                
                # plotting model slices 
                else:
                    ax.imshow(ims[k], aspect = 'auto', cmap = "jet", vmin = vmin, vmax = vmax)    

                    ax.plot(xSlices[k][0], xSlices[k][1], xSlices[k][2], linewidth = 0.5)
                    ax.plot(ySlices[k][0], ySlices[k][1], ySlices[k][2], linewidth = 0.5)
                    
                    ax.scatter(xshot[k], yshot[k], s = 20.0, color = "brown")
                    ax.scatter(xnode[k], ynode[k], s = 20.0, color = "gray")

                    ax.tick_params(direction = xTickDirection[k], axis='x') 
                    ax.tick_params(direction = yTickDirection[k], axis='y') 
                    
                    ax.set_xticks(xTickLock[k])
                    ax.set_yticks(yTickLock[k])

                    ax.set_xticklabels(xTickLabel[k])
                    ax.set_yticklabels(yTickLabel[k])
 
                    ax.set_xlabel(xLabel[k])
                    ax.set_ylabel(yLabel[k])
                    
                    if yInvert[k]:
                       ax.invert_yaxis()
    
    return None

#-------------------------------------------------------------------------

nx = 881
ny = 881
nz = 201

dx = 25
dy = 25
dz = 25

model = readBinaryVolume(nz, nx, ny, f"../inputs/models/testModel_{nz}x{nx}x{ny}_{dx}m.bin")

shots_file = "../inputs/geometry/xyz_shot_positions.txt"
nodes_file = "../inputs/geometry/xyz_node_positions.txt"

shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([nz/2, nx/2, ny/2], dtype = int) # [xy, zy, zx]
dh = np.array([dx, dy, dz])

check_geometry(model, shots, nodes, dh, slices, subplots, 2.8)
plt.savefig(f"modelTest.png", dpi = 200)

#-------------------------------------------------------------------------

v = np.array([1500, 2000, 3000, 4000])
z = np.array([1000, 1500, 2000])

x = np.sqrt((nodes[:,0] - shots[0])**2 + (nodes[:,1] - shots[1])**2)

fba = analytical_firstArrival(v, z, x)

n = 1256

dh = np.array([100, 50, 25], dtype = int)

pod = np.zeros((len(dh), n))
fim = np.zeros((len(dh), n))
fsm = np.zeros((len(dh), n))

for i in range(len(dh)):
    
    pod[i] = np.fromfile(f"../outputs/first_arrivals/{dh[i]}m_pod_data_nRec1256_shot_1.bin", dtype = np.float32, count = n)
    fim[i] = np.fromfile(f"../outputs/first_arrivals/{dh[i]}m_fim_data_nRec1256_shot_1.bin", dtype = np.float32, count = n)
    fsm[i] = np.fromfile(f"../outputs/first_arrivals/{dh[i]}m_fsm_data_nRec1256_shot_1.bin", dtype = np.float32, count = n)

offset = np.arange(n)

colors = ["blue", "orange", "green"]
styles = ["dashed", "dotted", "solid"]
titles = ["Podvin & Lecomte (1991)", "Jeong & Whitaker (2008)", "Detrixhe et al. (2013) | Noble et al. (2014)"]

xloc = np.linspace(0, n, 11, dtype = int)

fig, ax = plt.subplots(nrows = 3, ncols = 2, figsize = (15,8))

ax[0,0].plot(fba, color = "black")
ax[1,0].plot(fba, color = "black")
ax[2,0].plot(fba, color = "black")

for k in range(len(dh)):
    ax[0,0].plot(offset, pod[k], linestyle = styles[k], color = colors[0])
    ax[1,0].plot(offset, fim[k], linestyle = styles[k], color = colors[1])
    ax[2,0].plot(offset, fsm[k], linestyle = styles[k], color = colors[2])
    ax[0,1].plot(offset, fba - pod[k], linestyle = styles[k], color = colors[0])
    ax[1,1].plot(offset, fba - fim[k], linestyle = styles[k], color = colors[1])
    ax[2,1].plot(offset, fba - fsm[k], linestyle = styles[k], color = colors[2])

    for i in range(len(dh)):
        ax[i,0].set_xlabel("Trace index", fontsize= 15)
        ax[i,0].set_ylabel("Time [s]", fontsize = 15)
        ax[i,0].set_title(titles[i], fontsize = 18)

        ax[i,1].set_xlabel("Trace index", fontsize = 15)
        ax[i,1].set_ylabel("Diff = Ta - Tn [s]", fontsize = 15)
        ax[i,1].set_title(titles[i], fontsize = 18)

    for i in range(2):
        ax[k,i].set_xticks(xloc)
        ax[k,i].set_xticklabels(xloc)
        ax[k,i].set_xlim([0,n])

    ax[k,0].set_ylim([5.6, 6.4])
    ax[k,0].invert_yaxis()

plt.tight_layout()
plt.savefig(f"accuracyTest.png", dpi = 200)

#-------------------------------------------------------------------------

benchmark = np.loadtxt("runTime.txt", delimiter = ";", comments = "#")

bench_pod = benchmark[:3]
bench_fim = benchmark[3:6]
bench_fsm = benchmark[6:]

yaxis = ["Elapsed time [s]", "RAM usage [MB]", "GPU memory usage [MB]"]

xloc = [0, 1, 2]
xlab = ["2.9", "19.6", "156.1"]

fig, ax = plt.subplots(nrows = 3,ncols = 1, figsize = (10,9))

for i in range(len(dh)):
    ax[i].plot(bench_pod[:,i], "o--", label = "Podvin & Lecomte (1991)")
    ax[i].plot(bench_fim[:,i], "o--", label = "Jeong & Whitaker (2008)")
    ax[i].plot(bench_fsm[:,i], "o--", label = "Detrixhe et al. (2013) | Noble et al. (2014)")

    ax[i].set_xticks(xloc)
    ax[i].set_xticklabels(xlab)

    ax[i].legend(loc = "upper left")
    ax[i].set_ylabel(yaxis[i], fontsize = 15)
    ax[i].set_xlabel("Total samples in model [x 10⁶]", fontsize = 15)

plt.tight_layout()
plt.savefig(f"benchmarkTest.png", dpi = 200)
