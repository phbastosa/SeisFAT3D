import numpy as np
import segyio as sgy

from sys import argv

def catch_parameter(filename, target):
    file = open(filename,'r')
    for line in file.readlines():
        if line[0] != '#':
            splitted = line.split()
            if len(splitted) != 0:
                if splitted[0] == target: 
                    return splitted[2]

def read_binary_matrix(dim1, dim2, filename):
    data = np.fromfile(filename, dtype = np.float32, count = dim1*dim2)
    return np.reshape(data, [dim1, dim2], order = "F")

def show_binary_header(data):
    binHeader = sgy.binfield.keys
    print("\n Checking binary header \n")
    print(f"{'key': >25s} {'byte': ^6s} {'value': ^7s} \n")
    for k, v in binHeader.items():
        if v in data.bin:
            print(f"{k: >25s} {str(v): ^6s} {str(data.bin[v]): ^7s} \n")

def show_trace_header(data):
    traceHeader = sgy.tracefield.keys
    print("\n Checking trace header \n")
    print(f"{'Trace header': >40s} {'byte': ^6s} {'first': ^11s} {'last': ^11s} \n")
    for k, v in traceHeader.items():
        if v in data.bin:
            first = data.attributes(v)[0][0]
            last = data.attributes(v)[data.tracecount-1][0]
            print(f"{k: >40s} {str(v): ^6s} {str(first): ^11s} {str(last): ^11s}\n")

#---------------------------------------------------------------------

file = argv[1]

nt = int(catch_parameter(file, "time_samples"))
dt = float(catch_parameter(file, "time_spacing"))

shots_path = catch_parameter(file, "shots_file")
nodes_path = catch_parameter(file, "nodes_file")

shots = np.loadtxt(shots_path, dtype = float, delimiter = ",")
nodes = np.loadtxt(nodes_path, dtype = float, delimiter = ",")

total_shots = len(shots)
total_nodes = len(nodes)

for gather in range(total_nodes):

    offset = np.sqrt((nodes[gather, 0] - shots[:, 0])**2 + (nodes[gather, 1] - shots[:, 1])**2 + (nodes[gather, 2] - shots[:, 2])**2)  

    cmpx = 0.5 * (nodes[gather, 0] + shots[:, 0]) 
    cmpy = 0.5 * (nodes[gather, 1] + shots[:, 1]) 
    cmpt = np.sqrt(cmpx*cmpx + cmpy*cmpy) 

    fldr = 1001 + gather*np.ones(total_shots, dtype = int)

    data_path = f"../inputs/data/synthetic_seismogram_{nt}x{total_shots}_shot_{gather+1}.bin"
    segy_path = f"../inputs/data/seismogram_gather_{gather+1}.segy"

    data = read_binary_matrix(nt, total_shots, data_path)

    sgy.tools.from_array2D(segy_path, data.T)

    segy = sgy.open(segy_path, "r+", ignore_geometry = True)

    segy.bin[sgy.BinField.Interval]              = int(dt*1e6)
    segy.bin[sgy.BinField.IntervalOriginal]      = int(dt*1e6)
    segy.bin[sgy.BinField.Format]                = 1
    segy.bin[sgy.BinField.SortingCode]           = 1
    segy.bin[sgy.BinField.MeasurementSystem]     = 1
    segy.bin[sgy.BinField.ImpulseSignalPolarity] = 1

    tracl = 1 + np.arange(total_shots, dtype = int) 

    print(f"Convering {gather+1} of {total_nodes} .bin to .segy files...")

    for idx, key in enumerate(segy.header):

        key.update({sgy.TraceField.TRACE_SEQUENCE_LINE     : int(tracl[idx])          })
        key.update({sgy.TraceField.TRACE_SEQUENCE_FILE     : int(fldr[idx])           })
        key.update({sgy.TraceField.TRACE_SAMPLE_INTERVAL   : int(dt*1e6)              })
        key.update({sgy.TraceField.TRACE_SAMPLE_COUNT      : nt                       })
        key.update({sgy.TraceField.FieldRecord             : int(fldr[idx])           })
        key.update({sgy.TraceField.TraceNumber             : int(tracl[idx])          })
        key.update({sgy.TraceField.CDP                     : int(cmpt[idx]*1e2)       })
        key.update({sgy.TraceField.CDP_TRACE               : int(cmpt[idx]*1e2)       })
        key.update({sgy.TraceField.INLINE_3D               : int(shots[idx,0]*1e2)    })
        key.update({sgy.TraceField.CROSSLINE_3D            : int(shots[idx,1]*1e2)    })
        key.update({sgy.TraceField.TraceIdentificationCode : int(tracl[idx])          })
        key.update({sgy.TraceField.offset                  : int(offset[idx]*1e2)     })
        key.update({sgy.TraceField.ReceiverGroupElevation  : int(shots[idx,2]*1e2)    })
        key.update({sgy.TraceField.SourceSurfaceElevation  : int(nodes[gather,2]*1e2) })
        key.update({sgy.TraceField.ElevationScalar         : 100                      })
        key.update({sgy.TraceField.SourceGroupScalar       : 100                      })
        key.update({sgy.TraceField.SourceX                 : int(nodes[gather,0]*1e2) })
        key.update({sgy.TraceField.SourceY                 : int(nodes[gather,1]*1e2) })
        key.update({sgy.TraceField.SourceDepth             : int(nodes[gather,2]*1e2) })
        key.update({sgy.TraceField.GroupX                  : int(shots[idx,0]*1e2)    })
        key.update({sgy.TraceField.GroupY                  : int(shots[idx,1]*1e2)    })
        key.update({sgy.TraceField.GroupWaterDepth         : int(shots[idx,2]*1e2)    })
        key.update({sgy.TraceField.CoordinateUnits         : 1                        })
        key.update({sgy.TraceField.GainType                : 1                        })
        key.update({sgy.TraceField.TimeBaseCode            : 1                        })
        key.update({sgy.TraceField.CDP_X                   : int(cmpx[idx]*1e2)       })
        key.update({sgy.TraceField.CDP_Y                   : int(cmpy[idx]*1e2)       })
    
    segy.close()