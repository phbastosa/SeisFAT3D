import sys
import numpy as np
import matplotlib.pyplot as plt

n_inline = 51
n_crossline = 11

dc = 25
di = 12.5

offset_min = 25

nx_shots = 4
ny_shots = 4

x = np.linspace(2000, 8000, nx_shots)
y = np.linspace(2000, 8000, ny_shots)

total_shots = nx_shots * ny_shots
total_spreads = n_crossline * n_inline
total_stations = total_shots * total_spreads

x_spread = np.zeros(total_spreads)
y_spread = np.zeros(total_spreads)

shots = np.zeros((total_shots, 3))
nodes = np.zeros((total_stations, 3))

count = 0
for orientation in ["S-N", "N-S", "E-W", "W-E"]:
    
    count += 1
    for shot_index in range(total_shots):

        i = int(shot_index % ny_shots)
        j = int(shot_index / ny_shots)

        shots[i + j*ny_shots, 0] = x[i]
        shots[i + j*ny_shots, 1] = y[j] 

        for inline in range(n_inline):

            fill = slice(inline*n_crossline, inline*n_crossline + n_crossline)

            if orientation == "N-S":                
                y_spread[fill] = y[j] + offset_min + inline*di
                x_spread[fill] = np.linspace(x[i] - 0.5*n_crossline*dc, x[i] + 0.5*n_crossline*dc, n_crossline)

            elif orientation == "S-N":        
                y_spread[fill] = y[j] - offset_min - inline*di
                x_spread[fill] = np.linspace(x[i] - 0.5*n_crossline*dc, x[i] + 0.5*n_crossline*dc, n_crossline)
            
            elif orientation == "W-E": 
                x_spread[fill] = x[i] - offset_min - inline*di
                y_spread[fill] = np.linspace(y[j] - 0.5*n_crossline*dc, y[j] + 0.5*n_crossline*dc, n_crossline)

            elif orientation == "E-W":
                x_spread[fill] = x[i] + offset_min + inline*di
                y_spread[fill] = np.linspace(y[j] - 0.5*n_crossline*dc, y[j] + 0.5*n_crossline*dc, n_crossline)

            else:
                print("Orientation undefined!")
                print("Possibilities: N-S, S-N, W-E and E-W")
                sys.exit() 

        nodes[shot_index*total_spreads:shot_index*total_spreads + total_spreads, 0] = x_spread
        nodes[shot_index*total_spreads:shot_index*total_spreads + total_spreads, 1] = y_spread





    plt.figure(count, figsize = (6,6))
    plt.plot(nodes[:,0], nodes[:,1], "o", markersize = 1)
    plt.plot(shots[:,0], shots[:,1], "o")

    plt.xlim([0,10000])
    plt.ylim([0,10000])

plt.show()