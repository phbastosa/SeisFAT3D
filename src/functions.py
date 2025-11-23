import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def catch_parameter(parameters, target):
    file = open(parameters, "r")
    for line in file.readlines():
        if line[0] != "#":
            splitted = line.split()
            if len(splitted) != 0:
                if splitted[0] == target: 
                    return splitted[2]     

def read_binary_array(n1,filename):
    return np.fromfile(filename, dtype = np.float32, count = n1)    

def read_binary_matrix(n1,n2,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2)   
    return np.reshape(data, [n1, n2], order='F')

def read_binary_volume(n1, n2, n3, filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1, n2, n3], order = 'F')

def get_analytical_refractions(v, z, x):

    refracted_waves = np.zeros((len(z), len(x)))

    for n in range(len(z)):
        refracted_waves[n,:] = x / v[n+1]
        for i in range(n+1):
            angle = np.arcsin(v[i] / v[n+1])
            refracted_waves[n,:] += 2.0*z[i]*np.cos(angle) / v[i]
    
    return refracted_waves

def get_analytical_time_reflections(v, z, x):

    Tint = 2.0 * z * v[:-1]
    Vrms = np.zeros(len(z))

    reflections = np.zeros((len(z), len(x)))
    for i in range(len(z)):
        Vrms[i] = np.sqrt(np.sum(v[:i+1]**2 * Tint[:i+1]) / np.sum(Tint[:i+1]))
        reflections[i] = np.sqrt(x**2 + 4.0*np.sum(z[:i+1])**2) / Vrms[i]

    return reflections 

def get_analytical_amps_reflections(vp, ro, z):
    
    reflectivity = np.zeros(len(z))    
    reflectivity = (vp[1:]*ro[1:] - vp[:-1]*ro[:-1]) / (vp[1:]*ro[1:] + vp[:-1]*ro[:-1])

    reflections = np.zeros_like(reflectivity)
    transmission = np.zeros_like(reflectivity)

    reflections[0] = reflectivity[0]
    transmission[0] = 1.0 - reflectivity[0]

    for i in range(1, len(reflectivity)):
        reflections[i] = transmission[i-1]*reflectivity[i]
        transmission[i] = transmission[i-1]*(1.0 - reflectivity[i])

        for j in range(i,0,-1):
            reflections[i] *= 1.0 - reflectivity[i - j]

    return reflections

def get_ricker_wavelet(nt, dt, fmax):
    fc = fmax / (3.0*np.sqrt(np.pi))
    arg = np.pi*((np.arange(nt) - 0.5*nt)*dt*fc*np.pi)**2
    return (1.0 - 2.0*arg)*np.exp(-arg)

def compute_stiffness(vp, vs, ro, ep, dl, tht):
    
    SI = 1e9

    C = np.zeros((3,3))
    M = np.zeros((3,3))

    c11 = 0; c13 = 0; c15 = 0
    c33 = 0; c35 = 0; c55 = 0

    SI = 1e9

    c33 = ro*vp**2 / SI
    c55 = ro*vs**2 / SI

    c11 = c33*(1.0 + 2.0*ep)

    c13 = np.sqrt((c33 - c55)**2 + 2.0*dl*c33*(c33 - c55)) - c55

    C[0,0] = c11; C[0,1] = c13; C[0,2] = c15  
    C[1,0] = c13; C[1,1] = c33; C[1,2] = c35  
    C[2,0] = c15; C[2,1] = c35; C[2,2] = c55     

    tht = np.radians(tht)

    c = np.cos(tht)
    s = np.sin(tht)

    sin2 = np.sin(2.0*tht)
    cos2 = np.cos(2.0*tht)

    M = np.array([[     c**2,     s**2, sin2],
                  [     s**2,     c**2,-sin2],
                  [-0.5*sin2, 0.5*sin2, cos2]])
    
    Cr = (M @ C @ M.T) * SI

    return Cr

def plot_model_3D(model, dh, slices, **kwargs):

    m2km = 1e-3
    dbar = kwargs.get("dbar") if "dbar" in kwargs else 1.7
    scale = kwargs.get("scale") if "scale" in kwargs else 2.8 
    cmap = kwargs.get("cmap") if "cmap" in kwargs else "jet" 
    adjx = kwargs.get("adjx") if "adjx" in kwargs else 0.5
    cblab = kwargs.get("cblab") if "cblab" in kwargs else "Velocity [km/s]"

    shots_defined = True if "shots" in kwargs else False
    nodes_defined = True if "nodes" in kwargs else False
    eikonal_defined = True if "eikonal" in kwargs else False

    shots = np.loadtxt(kwargs.get("shots"), delimiter = ',') if shots_defined else None
    nodes = np.loadtxt(kwargs.get("nodes"), delimiter = ',') if nodes_defined else None

    if eikonal_defined:
        eikonal = kwargs.get("eikonal")
        eikonal_levels = kwargs.get("eikonal_levels")
        eikonal_colors = kwargs.get("eikonal_colors")

    model_shape = np.array(np.shape(model))

    max_model_distance = np.max(model_shape)
    min_model_distance = np.min(model_shape)

    vmin = np.min(model) if kwargs.get("vmin") == None else kwargs.get("vmin") 
    vmax = np.max(model) if kwargs.get("vmax") == None else kwargs.get("vmax") 
    
    nz, nx, ny = model_shape

    [z, x, y] = scale*(min_model_distance / max_model_distance) * model_shape / max_model_distance

    px = 1.0 / plt.rcParams['figure.dpi']  
    ticks = np.array([3, 7, 7], dtype = int)

    fig = plt.figure(figsize = (900*px, 780*px))

    xloc = np.linspace(0, nx-1, ticks[1], dtype = int)
    yloc = np.linspace(0, ny-1, ticks[2], dtype = int)
    zloc = np.linspace(0, nz-1, ticks[0], dtype = int)

    xlab = np.around(xloc * dh[0] * m2km, decimals = 1)
    ylab = np.around(yloc * dh[1] * m2km, decimals = 1)
    zlab = np.around(zloc * dh[2] * m2km, decimals = 1)

    axes = np.array([[adjx - x, 0.98 - y      , x, y], 
                     [    adjx, 0.98 - y      , z, y],
                     [adjx - x, 0.98 - y - z  , x, z],
                     [adjx - x, 0.98 - y - dbar*z, x, z]])

    xTickDirection = ['out', 'out', 'out']
    yTickDirection = ['out', 'in', 'out']

    xTickLock = [xloc, zloc[1:], xloc]
    yTickLock = [yloc, yloc, zloc[1:]]

    xTickLabel = [[], zlab[1:], xlab]
    yTickLabel = [ylab, [], zlab[1:]]

    xLabel = ["X [km]", "Z [km]", "X [km]"]
    yLabel = ["Y [km]", "      ", "Z [km]"]

    yInvert = [True, True, False]

    xSlices = [[np.arange(model_shape[1]), np.ones(model_shape[1])*slices[1], "--g"],
               [np.arange(model_shape[0]), np.ones(model_shape[0])*slices[1], "--g"],
               [np.arange(model_shape[1]), np.ones(model_shape[1])*slices[0], "--r"]] 

    ySlices = [[np.ones(model_shape[2])*slices[2], np.arange(model_shape[2]), "--m"],
               [np.ones(model_shape[2])*slices[0], np.arange(model_shape[2]), "--r"],
               [np.ones(model_shape[0])*slices[2], np.arange(model_shape[0]), "--m"]]

    if shots_defined:
    
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
            
    if nodes_defined:

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

    if eikonal_defined:    

        eiks = []

        if len(np.shape(eikonal)) == 4:

            for i in range(len(eikonal)):
                curve = [eikonal[i,slices[0],:,:].T, eikonal[i,:,slices[2],:].T, eikonal[i,:,:,slices[1]]] 
                eiks.append(curve)
        else: 
            
            eiks = [eikonal[slices[0],:,:].T, eikonal[:,slices[2],:].T, eikonal[:,:,slices[1]]]

    ims = [model[slices[0],:,:].T, model[:,slices[2],:].T, model[:,:,slices[1]]]

    if shots_defined:

        if np.size(shots) != 3:
            xshot = [shots[:,0]/dh[0], zy_plane_shot_z, zx_plane_shot_x]
            yshot = [shots[:,1]/dh[1], zy_plane_shot_y, zx_plane_shot_z]
        else:
            xshot = [shots[0]/dh[0], zy_plane_shot_z, zx_plane_shot_x]
            yshot = [shots[1]/dh[1], zy_plane_shot_y, zx_plane_shot_z]

    if nodes_defined:

        if np.size(nodes) != 3:
            xnode = [nodes[:,0]/dh[0], zy_plane_node_z, zx_plane_node_x]
            ynode = [nodes[:,1]/dh[1], zy_plane_node_y, zx_plane_node_z]
        else:
            xnode = [nodes[0]/dh[0], zy_plane_node_z, zx_plane_node_x]
            ynode = [nodes[1]/dh[1], zy_plane_node_y, zx_plane_node_z]

    for k, axs in enumerate(axes):      
        
        ax = fig.add_axes(axs)                         
        
        if k == 3:

            ax.axis("off")
            cmap = mpl.colormaps[cmap]
            norm = mpl.colors.Normalize(vmin*m2km, vmax*m2km)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("bottom", size="5%", pad=0)
            cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap), cax = cax, ticks = np.linspace(vmin*m2km, vmax*m2km, 5), orientation = "horizontal")
            cbar.ax.set_xticklabels(np.around(np.linspace(vmin*m2km, vmax*m2km, 5), decimals = 1))
            cbar.set_label(cblab, fontsize = 15)
         
        else:
            
            ax.imshow(ims[k], aspect = 'auto', cmap = cmap, vmin = vmin, vmax = vmax)    
            
            ax.plot(xSlices[k][0], xSlices[k][1], xSlices[k][2], linewidth = 0.5)
            ax.plot(ySlices[k][0], ySlices[k][1], ySlices[k][2], linewidth = 0.5)
            
            if eikonal_defined:
                ax.contour(eiks[k], levels = eikonal_levels, colors = eikonal_colors, linestyles = "dashed")

            if nodes_defined:
                ax.plot(xnode[k], ynode[k], "o", markersize = 1, color = "gray")
            
            if shots_defined:
                ax.plot(xshot[k], yshot[k], "*", markersize = 5, color = "green")

            ax.tick_params(direction = xTickDirection[k], axis='x') 
            ax.tick_params(direction = yTickDirection[k], axis='y') 
            
            ax.set_xticks(xTickLock[k])
            ax.set_yticks(yTickLock[k])

            ax.set_xticklabels(xTickLabel[k])
            ax.set_yticklabels(yTickLabel[k])

            ax.set_xlabel(xLabel[k], fontsize = 15)
            ax.set_ylabel(yLabel[k], fontsize = 15)
            
            if yInvert[k]:
                ax.invert_yaxis()
    
    return None    
