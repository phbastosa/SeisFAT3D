import numpy as np

from scipy.signal import savgol_filter

nline_shots = 10
nline_nodes = 51

total_shots = nline_shots*nline_shots 
total_nodes = nline_nodes*nline_nodes

dt = 1e-3
fmax = 25
raw_nt = 6001
time_cut = 5.0

t0 = int(np.pi / fmax / dt) 
tmax = int(time_cut / dt) + 1

filter_window = 9

for shot in range(total_shots):
    
    raw_picks = np.fromfile(f"../inputs/data/rawPicks_{total_nodes}_shot_{shot+1}.bin", dtype = np.float32, count = total_nodes)

    raw_picks = raw_picks - 0.5*t0*dt

    picks = np.zeros(total_nodes)

    for line in range(nline_nodes):

        window = slice(int(line*total_nodes/nline_nodes), int((line+1)*total_nodes/nline_nodes))
    
        aux_picks = raw_picks[window].copy()

        picks[window] = savgol_filter(aux_picks, filter_window, 2) 

    print(f"File ../inputs/data/obsData_nRec{total_nodes}_shot_{shot+1}.bin was written successfully.")
    picks.astype("float32", order = "F").tofile(f"../inputs/data/obsData_nRec{total_nodes}_shot_{shot+1}.bin")


