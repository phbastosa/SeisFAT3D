import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def readBinaryVolume(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

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

#---------------------------------------------------------------------

nx = 101
ny = 101
nz = 41

dx = 50
dy = 50
dz = 50

model = readBinaryVolume(nz, nx, ny, f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dx}m.bin")

shots_file = "../inputs/geometry/xyz_shot_positions.txt"
nodes_file = "../inputs/geometry/xyz_node_positions.txt"

shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([nz/2, nx/2, ny/2], dtype = int) # [xy, zy, zx]
dh = np.array([dx, dy, dz])

check_geometry(model, shots, nodes, dh, slices, subplots, 1.2)
plt.savefig(f"trueModel.png", dpi = 200)

check_geometry(1500.0*np.ones_like(model), shots, nodes, dh, slices, subplots, 1.2)
plt.savefig(f"initModel.png", dpi = 200)
# --------------------------------------------------

convergence_ls = np.loadtxt("../outputs/convergence/ls_convergence_3_iterations.txt")
convergence_adj = np.loadtxt("../outputs/convergence/adj_convergence_3_iterations.txt")

plt.figure(2, figsize = (10,4))
plt.plot(convergence_ls, "o--", label = "Least squares approach")
plt.plot(convergence_adj, "o--", label = "Adjoint state approach")

plt.title("Convergence curve", fontsize = 18)
plt.xlabel("Iteration number", fontsize = 15)
plt.ylabel("Objective function L2 norm", fontsize = 15)

plt.legend(loc = "upper right", fontsize = 12)
plt.tight_layout()
plt.savefig(f"curve.png", dpi = 200)

#-----------------------------------------------

final_model_ls = readBinaryVolume(nz, nx, ny, f"../outputs/recovered_models/ls_final_model_{nz}x{nx}x{ny}.bin")
final_model_adj = readBinaryVolume(nz, nx, ny, f"../outputs/recovered_models/adj_final_model_{nz}x{nx}x{ny}.bin")

check_geometry(final_model_ls, shots, nodes, dh, slices, subplots, 1.2)
plt.savefig(f"final_model_ls.png", dpi = 200)

check_geometry(final_model_adj, shots, nodes, dh, slices, subplots, 1.2)
plt.savefig(f"final_model_adj.png", dpi = 200)

