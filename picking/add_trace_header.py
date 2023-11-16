import numpy as np
import segyio as sgy
import matplotlib.pyplot as plt

def read_binary_matrix(dim1, dim2, filename):
    data = np.fromfile(filename, dtype = np.float32, count = dim1*dim2)
    return np.reshape(data, [dim1, dim2], order = "F")

shots_path = "../inputs/geometry/xyz_shot_positions.txt"
nodes_path = "../inputs/geometry/xyz_node_positions.txt"

shots = np.loadtxt(shots_path, dtype = float, delimiter = ",")
nodes = np.loadtxt(nodes_path, dtype = float, delimiter = ",")

total_shots = len(shots)
total_nodes = len(nodes)

offset = np.sqrt((shots[0,0] - nodes[:,0])**2 + (shots[0,1] - nodes[:,1])**2 + (shots[0,2] - nodes[:,2])**2)  

plt.plot(shots[:,0], shots[:,1], "o")
plt.plot(nodes[:,0], nodes[:,1], "o")
plt.show()