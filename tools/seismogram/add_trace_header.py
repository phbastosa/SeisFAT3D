import numpy as np
import segyio as sgy

def read_binary_matrix(dim1, dim2, filename):
    data = np.fromfile(filename, dtype = np.float32, count = dim1*dim2)
    return np.reshape(data, [dim1, dim2], order = "F")

nt = 5001

shots_path = "../../inputs/geometry/xyz_shot_positions.txt"
nodes_path = "../../inputs/geometry/xyz_node_positions.txt"

data_path = "data_bin/output_seismogram_5001x1681_shot_1.bin"
segy_path = "data_sgy/seismogram_shot_1.segy"

shots = np.loadtxt(shots_path, dtype = float, delimiter = ",")
nodes = np.loadtxt(nodes_path, dtype = float, delimiter = ",")

total_shots = len(shots)
total_nodes = len(nodes)

offset = np.sqrt((nodes[0,0] - shots[:,0])**2 + (nodes[0,1] - shots[:,1])**2 + (nodes[0,2] - shots[:,2])**2)  

cmpx = 0.5 * (nodes[0,0] + shots[:,0]) 
cmpy = 0.5 * (nodes[0,1] + shots[:,1]) 

data = read_binary_matrix(nt, total_shots, data_path)

sgy.tools.from_array2D(segy_path, data.T)

segy = sgy.open(segy_path, "r+", ignore_geometry = True)
