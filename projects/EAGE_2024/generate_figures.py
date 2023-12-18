import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def readBinaryVolume(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

def readBinaryMatrix(n1,n2,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2)    
    return np.reshape(data, [n1,n2], order='F')

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

def check_full_model(models, shots, nodes, dh, slices, subplots, vmin, vmax, scale):
    
    if np.sum(subplots) == 2:
        modelShape = np.array(np.shape(models))
        maxModelDistance = np.max(np.shape(models))
        minModelDistance = np.min(np.shape(models))

    else:
        modelShape = np.array(np.shape(models[0]))
        maxModelDistance = np.max(np.shape(models[0]))
        minModelDistance = np.min(np.shape(models[0]))

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
                     [0.75 - x, 0.98 - y - 1.5*z, x, z]])

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
                imPod = [volPod[slices[0],:,:].T, volPod[:,slices[2],:].T, volPod[:,:,slices[1]]]
                imFim = [volFIM[slices[0],:,:].T, volFIM[:,slices[2],:].T, volFIM[:,:,slices[1]]]
                imFsm = [volFSM[slices[0],:,:].T, volFSM[:,slices[2],:].T, volFSM[:,:,slices[1]]]           
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
                    
                    ax.contour(imPod[k], levels = [1.0, 2.0, 3.0], colors = "blue", linestyles = "dashed")    
                    ax.contour(imFim[k], levels = [1.0, 2.0, 3.0], colors = "orange", linestyles = "dashed")    
                    ax.contour(imFsm[k], levels = [1.0, 2.0, 3.0], colors = "green", linestyles = "dashed")    

                    ax.plot(xSlices[k][0], xSlices[k][1], xSlices[k][2], linewidth = 0.5)
                    ax.plot(ySlices[k][0], ySlices[k][1], ySlices[k][2], linewidth = 0.5)
                    
                    ax.scatter(xshot[k], yshot[k], s = 20.0, color = "black")
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

#---------------------------------------------------------------------

nx = 441
ny = 441
nz = 161

dh = 25

# accuracy_model = readBinaryVolume(nz, nx, ny, f"../inputs/models/accuracyModelTest_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

shots_file = "../inputs/geometry/xyz_shots_position.txt"
nodes_file = "../inputs/geometry/xyz_nodes_position.txt"
 
shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([int(nz/2), int(nx/2), 120], dtype = int) # [xy, zy, zx]
dh = np.array([dh, dh, dh])

vmin = 2000
vmax = 5000

# volPod = readBinaryVolume(nz, nx, ny, f"../outputs/travel_times/pod_time_volume_{nz}x{nx}x{ny}_shot_1.bin")
# volFIM = readBinaryVolume(nz, nx, ny, f"../outputs/travel_times/fim_time_volume_{nz}x{nx}x{ny}_shot_1.bin")
# volFSM = readBinaryVolume(nz, nx, ny, f"../outputs/travel_times/fsm_time_volume_{nz}x{nx}x{ny}_shot_1.bin")

# check_full_model(accuracy_model, shots, nodes, dh, slices, subplots, vmin, vmax, 1.5)
# plt.savefig(f"accuracy_model_test.png", dpi = 200)
# plt.clf()

nt = 6001
dt = 1e-3

nTraces = len(nodes)

velocity = np.array([2000, 3000, 4000, 5000])
thickness = np.array([1000, 1000, 1000])

offset = np.sqrt((shots[0] - nodes[:,0])**2 + (shots[1] - nodes[:,1])**2 + (shots[2] - nodes[:,2])**2)

analytical_times = analytical_firstArrival(velocity, thickness, offset)

pod_firstArrivals = np.fromfile(f"../outputs/first_arrivals/pod_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)
fim_firstArrivals = np.fromfile(f"../outputs/first_arrivals/fim_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)
fsm_firstArrivals = np.fromfile(f"../outputs/first_arrivals/fsm_data_nRec{nTraces}_shot_1.bin", dtype = np.float32, count = nTraces)

t0 = int(np.pi / 25.0 / dt) 
tmax = int(5.0 / dt) + 1

seismic = readBinaryMatrix(nt, nTraces, f"../inputs/data/synthetic_seismogram_6001x273_shot_1.bin")

seismic = seismic[t0:tmax+t0,:]

sl1 = slice(0, int(nTraces/3))
sl2 = slice(int(nTraces/3), int(2*nTraces/3))
sl3 = slice(int(2*nTraces/3), nTraces)

slices = [sl1, sl2, sl3]

scale = 0.1*np.std(seismic)

tloc = np.linspace(0, tmax-1, 11)
tlab = np.around(tloc * dt, decimals = 3)

xloc = np.linspace(0, nTraces/3-1, 7)
xlab = np.linspace(0, nTraces/3, 7, dtype = int)

fig, ax = plt.subplots(ncols = 3, nrows = 1, figsize = (10,6))

for i in range(len(slices)):

    ax[i].imshow(seismic[:, slices[i]], aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
    ax[i].plot(analytical_times[slices[i]]/dt, "o", color = "red", markersize = 1, label = "Analytical travel times")
    ax[i].plot(pod_firstArrivals[slices[i]]/dt, "o", color = "blue", markersize = 1, label = "Podvin & Lecomte (1991)")
    ax[i].plot(fim_firstArrivals[slices[i]]/dt, "o", color = "orange", markersize = 1, label = "Jeong & Whitaker (2008)")
    ax[i].plot(fsm_firstArrivals[slices[i]]/dt, "o", color = "green", markersize = 1, label = "Noble, Gesret & Belayouni (2014)")

    ax[i].set_yticks(tloc)
    ax[i].set_yticklabels(tlab)

    ax[i].set_xticks(xloc)
    ax[i].set_xticklabels(xlab)

    ax[i].set_xlabel("Trace index", fontsize = 12)
    ax[i].set_ylabel("Time [s]", fontsize = 12)

    ax[i].legend(loc = "upper right", fontsize = 8)

fig.tight_layout()
plt.savefig("modeled_seismograms.png", dpi = 300)
plt.show()

