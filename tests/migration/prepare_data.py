import numpy as np
import matplotlib.pyplot as plt

from sys import path
path.append("../src/")
import functions

ns = 4
nr = 81
nt = 2001

dt = 1e-3

fmax = 30.0

tlag = 2.0*np.sqrt(np.pi) / fmax

input_data = np.zeros((nt, nr))

for i in range(ns):

    h_data_path = f"../outputs/seismograms/h_acoustic_seismogram_Nsamples2001_nRec81_shot_{i+1}.bin"
    d_data_path = f"../outputs/seismograms/d_acoustic_seismogram_Nsamples2001_nRec81_shot_{i+1}.bin"

    h_data = functions.read_binary_matrix(nt, nr, h_data_path)
    d_data = functions.read_binary_matrix(nt, nr, d_data_path)

    diff = d_data - h_data

    input_data[:-int(tlag/dt)] = diff[int(tlag/dt):]

    input_data.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/data/migration_test_data_nSamples{nt}_nRec{nr}_shot_{i+1}.bin")



