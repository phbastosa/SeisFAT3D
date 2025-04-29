import numpy as np
import segyio as sgy
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from matplotlib.widgets import Cursor
from mpl_toolkits.axes_grid1 import make_axes_locatable

def read_binary_array(n1, filename):
    return np.fromfile(filename, dtype = np.float32, count = n1)    

def read_binary_matrix(n1, n2, filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2)    
    return np.reshape(data, [n1, n2], order = 'F')

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
                ax.plot(xnode[k], ynode[k], "o", markersize = 5, color = "gray")
            
            if shots_defined:
                ax.plot(xshot[k], yshot[k], "*", markersize = 5, color = "black")

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


class Pick():
    def __init__(self):
        self.x = np.array([])
        self.t = np.array([])

class ManualPicking(Pick):
    
    def __init__(self, data_path, gain):
        
        self.gain = gain

        self.data = sgy.open(data_path, ignore_geometry = True)

        self.traces = [self.data.attributes(13)[self.data.tracecount-1][0]]

        self.nt = self.data.attributes(115)[0][0]
        self.dt = self.data.attributes(117)[0][0] * 1e-6

        self.sx = self.data.attributes(73)[:] * 1e-5
        self.sy = self.data.attributes(77)[:] * 1e-5

        self.gx = self.data.attributes(81)[:] * 1e-5
        self.gy = self.data.attributes(85)[:] * 1e-5

        self.seismic = self.data.trace.raw[:].T

        self.picks = []
        self.gathers = []

        for i in range(len(self.traces)):
            self.picks.append(Pick())
            self.gathers.append(self.seismic[:,sum(self.traces[:i]):sum(self.traces[:i+1])])

        self.current_plot = 0

        for i in range(len(self.traces)):
            while True:
                if self.current_plot == i:
                    print("Close plot windows to interpolate picks")
                    self.plot_figures()
                    self.interpolation()
                    self.plot_figures()

                    self.entry = str(input("Continue picking? (y or n) "))
                    if self.entry in "nN":
                        self.save_picks()
                        break

            self.current_plot += 1

    def plot_geometry(self):
        
        fig, ax = plt.subplots(dpi = 100, tight_layout = True)        
        manager = plt.get_current_fig_manager()
        manager.window.setGeometry(0,0,500,500) 

        ax.plot(self.gx, self.gy, "ob", markersize = 1, alpha = 0.6, label = "Nodes position")

        for i in range(len(self.traces)):
            if self.current_plot == i:
                trace_slicer = slice(sum(self.traces[:i]), sum(self.traces[:i+1]))
                ax.plot(self.gx[trace_slicer], self.gy[trace_slicer], "og", markersize = 3, label = "Active gather")

        ax.plot(self.sx, self.sy, "or", markersize = 5, label = "Shots position")

        ax.set_title("Acquisition Geometry")
        ax.set_xlabel("UTM E [km]")
        ax.set_ylabel("UTM N [km]")
        ax.legend(loc = "lower left")
        
        fig.tight_layout()

    def plot_seismic(self):

        self.fig, self.ax = plt.subplots(dpi = 100, tight_layout = True)

        tloc = np.linspace(0, self.nt, 21)
        tlab = np.around(tloc*self.dt, decimals = 2)

        manager = plt.get_current_fig_manager()
        manager.window.setGeometry(600,0,1200,900) 
        
        self.cursor = Cursor(self.ax, horizOn = True, vertOn = True, color = "green", linewidth = 1)

        self.ax.set_title(f"Gather {self.current_plot+1}")
        self.ax.set_xlabel("Shot index")
        self.ax.set_ylabel("Time [s]")

        for i in range(len(self.traces)):
            if self.current_plot == i:
                
                scale = self.gain * np.std(self.gathers[i])
                
                self.ax.imshow(self.gathers[i], aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
                self.ax.plot(self.picks[i].x, self.picks[i].t, "or")

                self.ax.set_yticks(tloc)
                self.ax.set_yticklabels(tlab)

                xloc = np.linspace(0, self.traces[i]-1, 9, dtype = int)
                xlab = np.linspace(1, self.traces[i], 9, dtype = int)

                self.ax.set_xticks(xloc)
                self.ax.set_xticklabels(xlab)

                self.fig.tight_layout()

                self.fig.canvas.mpl_connect("key_press_event", self.pick_data)

    def pick_data(self, event):

        for i in range(len(self.traces)):
            if self.current_plot == i:

                if event.key == "t":
                    self.picks[i].x = np.append(self.picks[i].x, np.abs(int(event.xdata)))
                    self.picks[i].t = np.append(self.picks[i].t, event.ydata)       

                    print(f"Pick (x,t) = ({self.picks[i].x[-1]:.0f}, {self.picks[i].t[-1]*self.dt:.3f}) picked")

                elif event.key == "d":
                    
                    print(f"Pick (x,t) = ({self.picks[i].x[-1]:.0f}, {self.picks[i].t[-1]*self.dt:.3f}) deleted")

                    self.picks[i].x = self.picks[i].x[:-1]
                    self.picks[i].t = self.picks[i].t[:-1]    
                    
                self.ax.plot(self.picks[i].x, self.picks[i].t, "or")
    
    def plot_figures(self):
        plt.close("all")
        self.plot_seismic()
        self.plot_geometry()
        plt.show()

    def interpolation(self):
        for i in range(len(self.traces)):
            
            if len(self.picks[i].x) < 3:
                return None
            
            if self.current_plot == i:
                
                x = np.arange(np.min(self.picks[i].x), np.max(self.picks[i].x) + 1)        

                f = interp1d(self.picks[i].x, self.picks[i].t, kind = "slinear")

                t = f(x)

                self.picks[i].x = x.copy()
                self.picks[i].t = t.copy()

    def save_picks(self):

        output_file = f"{self.picks_path}_shotGather_{self.current_plot+1}.txt"

        with open(output_file, "a") as file:
            for i in range(len(self.traces)):
                if self.current_plot == i:
                    for n in range(len(self.picks[i].x)):
                        file.write(f"{self.picks[i].x[n] + sum(self.traces[:i]):.0f}, {self.picks[i].t[n] * self.dt:.3f}\n")            
        
        file.close()
