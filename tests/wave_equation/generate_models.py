import numpy as np

from sys import path, argv
path.append("../src/")
import functions

nx = int(functions.catch_parameter(argv[1], "x_samples"))
ny = int(functions.catch_parameter(argv[1], "y_samples"))
nz = int(functions.catch_parameter(argv[1], "z_samples"))

dz = float(functions.catch_parameter(argv[1], "z_spacing"))

v = np.array([1500, 2000, 3000], dtype = float)
z = np.array((300, 400), dtype = float)

functions.build_layer_cake_model(z, v, nx, ny, nz, dz)

