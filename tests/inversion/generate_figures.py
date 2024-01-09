import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def readBinaryVolume(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

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

nx = 201
ny = 201
nz = 51

dx = 100
dy = 100
dz = 100

true_model = readBinaryVolume(nz, nx, ny, f"../inputs/models/trueModelTest_{nz}x{nx}x{ny}_{dx}m.bin")
init_model = readBinaryVolume(nz, nx, ny, f"../inputs/models/initModelTest_{nz}x{nx}x{ny}_{dx}m.bin")

shots_file = "../inputs/geometry/xyz_shots_position.txt"
nodes_file = "../inputs/geometry/xyz_nodes_position.txt"
 
shots = np.loadtxt(shots_file, delimiter = ',')
nodes = np.loadtxt(nodes_file, delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([20, 75, 75], dtype = int) # [xy, zy, zx]
dh = np.array([dx, dy, dz])

vmin = 2000
vmax = 4500

vdiff = 500

check_full_model(true_model, shots, nodes, dh, slices, subplots, vmin, vmax, 2.5)
plt.savefig(f"trueModel.png", dpi = 200)
plt.clf()

check_full_model(init_model, shots, nodes, dh, slices, subplots, vmin, vmax, 2.5)
plt.savefig(f"initModel.png", dpi = 200)
plt.clf()

check_full_model(true_model - init_model, shots, nodes, dh, slices, subplots, -vdiff, vdiff, 2.5)
plt.savefig(f"diffModel.png", dpi = 200)
plt.clf()

#-----------------------------------------------

final_model_ls = readBinaryVolume(nz, nx, ny, f"../outputs/recovered_models/ls_final_model_{nz}x{nx}x{ny}.bin")
final_model_adj = readBinaryVolume(nz, nx, ny, f"../outputs/recovered_models/adj_final_model_{nz}x{nx}x{ny}.bin")

check_full_model(final_model_ls, shots, nodes, dh, slices, subplots, vmin, vmax, 2.5)
plt.savefig(f"final_model_ls.png", dpi = 200)
plt.clf()

check_full_model(final_model_adj, shots, nodes, dh, slices, subplots, vmin, vmax, 2.5)
plt.savefig(f"final_model_adj.png", dpi = 200)
plt.clf()

check_full_model(final_model_ls - init_model, shots, nodes, dh, slices, subplots, -vdiff, vdiff, 2.5)
plt.savefig(f"diff_model_ls.png", dpi = 200)
plt.clf()

check_full_model(final_model_adj - init_model, shots, nodes, dh, slices, subplots, -vdiff, vdiff, 2.5)
plt.savefig(f"diff_model_adj.png", dpi = 200)
plt.clf()

# --------------------------------------------------

circles = np.array([[7500, 7500], [11500, 7500], [7500, 11500], [11500, 11500]]) / dx

logs = np.zeros((len(circles), 3, nz))

for i in range(len(circles)):
    logs[i,0,:] = true_model[:, int(circles[i,0]), int(circles[0,1])] - init_model[:, int(circles[i,0]), int(circles[0,1])]
    logs[i,1,:] = final_model_ls[:, int(circles[i,0]), int(circles[0,1])] - init_model[:, int(circles[i,0]), int(circles[0,1])]
    logs[i,2,:] = final_model_adj[:, int(circles[i,0]), int(circles[0,1])] - init_model[:, int(circles[i,0]), int(circles[0,1])]

depth = np.arange(nz)*dz

for i in range(len(circles)):

    plt.figure(i+1, figsize = (4,6))
    plt.plot(0.0*depth, depth, color = "black")
    plt.plot(logs[i,0,:], depth, color = "red")
    plt.plot(logs[i,1,:], depth, color = "orange")
    plt.plot(logs[i,2,:], depth, color = "green")

    plt.title(f"(x,y) = ({circles[i,0]*dx:.0f}, {circles[i,1]*dy:.0f}) m", fontsize = 15)
    plt.xlabel("Velocity anomaly [m/s]", fontsize = 15)
    plt.ylabel("Depth [m]", fontsize = 15)

    plt.xlim([-1000,1000])
    plt.ylim([0,5000])
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f"log{i+1}_test.png", dpi = 200)

plt.clf()

# --------------------------------------------------
# model correlation

ls_diff = final_model_ls - true_model
adj_diff = final_model_adj - true_model

rms_error_ls = np.sqrt(np.sum(ls_diff**2)/(nx*ny*nz))
rms_error_adj = np.sqrt(np.sum(adj_diff**2)/(nx*ny*nz))

max_error_ls = np.max(ls_diff)
max_error_adj = np.max(adj_diff)

min_error_ls = np.min(ls_diff)
min_error_adj = np.min(adj_diff)

print(f"|{'-'*74}|")
print(f"| {'Model difference analysis':^30s} | {'RMS ERROR':^11s} | {'MAX ERROR':^11s} | {'MIN ERROR':^11s} |")
print(f"|{'-'*74}|")
print(f"| {'Classical Tomography':^30s} | {f'{rms_error_ls:.4f}':^11s} | {f'{max_error_ls:.5f}':^11s} | {f'{min_error_ls:.4f}':^11s} |")
print(f"| {'Adjoint State Tomography':^30s} | {f'{rms_error_adj:.4f}':^11s} | {f'{max_error_adj:.5f}':^11s} | {f'{min_error_adj:.4f}':^11s} |")
print(f"|{'-'*74}|\n")

ns = len(shots)
nr = len(nodes)

ls_data = np.zeros(ns*nr)
adj_data = np.zeros(ns*nr)
obs_data = np.zeros(ns*nr)

for i in range(ns):
    ls_data[i*nr:nr+i*nr] = np.fromfile(f"../outputs/first_arrivals/ls_fsm_data_nRec400_shot_{i+1}.bin", count = nr, dtype = np.float32)
    adj_data[i*nr:nr+i*nr] = np.fromfile(f"../outputs/first_arrivals/adj_fsm_data_nRec400_shot_{i+1}.bin", count = nr, dtype = np.float32)
    obs_data[i*nr:nr+i*nr] = np.fromfile(f"../inputs/data/obs_fsm_data_nRec400_shot_{i+1}.bin", count = nr, dtype = np.float32)

rms_error_ls = np.sqrt(np.sum((obs_data - ls_data)**2)/(ns*nr))
rms_error_adj = np.sqrt(np.sum((obs_data - adj_data)**2)/(ns*nr))

max_error_ls = np.max(obs_data - ls_data)
max_error_adj = np.max(obs_data - adj_data)

min_error_ls = np.min(obs_data - ls_data)
min_error_adj = np.min(obs_data - adj_data)

print(f"|{'-'*74}|")
print(f"| {'Data difference analysis':^30s} | {'RMS ERROR':^11s} | {'MAX ERROR':^11s} | {'MIN ERROR':^11s} |")
print(f"|{'-'*74}|")
print(f"| {'Classical Tomography':^30s} | {f'{rms_error_ls:.4f}':^11s} | {f'{max_error_ls:.5f}':^11s} | {f'{min_error_ls:.4f}':^11s} |")
print(f"| {'Adjoint State Tomography':^30s} | {f'{rms_error_adj:.4f}':^11s} | {f'{max_error_adj:.5f}':^11s} | {f'{min_error_adj:.4f}':^11s} |")
print(f"|{'-'*74}|")

# --------------------------------------------------
convergence_ls = np.loadtxt("../outputs/convergence/ls_convergence_5_iterations.txt")
convergence_adj = np.loadtxt("../outputs/convergence/adj_convergence_5_iterations.txt")

plt.figure(21, figsize = (10,4))
plt.plot(convergence_ls, "o--", label = "Least squares approach", color = "orange")
plt.plot(convergence_adj, "o--", label = "Adjoint state approach", color = "green")

plt.title("Convergence curve", fontsize = 18)
plt.xlabel("Iteration number", fontsize = 15)
plt.ylabel("Objective function L2 norm", fontsize = 15)
 
plt.grid(True)
plt.tight_layout()
plt.savefig(f"curve.png", dpi = 200)
