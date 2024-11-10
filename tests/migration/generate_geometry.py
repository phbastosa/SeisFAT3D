import numpy as np

nsx = 41
nsy = 4

sx, sy = np.meshgrid(np.linspace(500, 2500, nsx), np.linspace(250, 1750, nsy))

nrx = 60 
nry = 40

rx, ry = np.meshgrid(np.linspace(25, 2975, nrx), np.linspace(25, 1975, nry))

spread_x = 20
spread_y = 10

SPS = np.zeros((nsx*nsy, 3), dtype = float)
RPS = np.zeros((nrx*nry, 3), dtype = float)
XPS = np.zeros((nsx*nsy, 3), dtype = int)

spread = spread_x*spread_y

rec_index = np.zeros((nry, nrx))

blk_count = 0

for k in range(0, nry, spread_y):
    for j in range(0, nrx, spread_x):

        xblk = slice(j, j + spread_x)    
        yblk = slice(k, k + spread_y)     
        
        blkId = slice(blk_count*spread, blk_count*spread + spread)

        RPS[blkId, 0] = np.reshape(rx[yblk, xblk], spread, order = "F")
        RPS[blkId, 1] = np.reshape(ry[yblk, xblk], spread, order = "F")

        rec_index[yblk, xblk] = np.reshape(np.arange(spread) + blk_count*spread, [spread_y, spread_x], order = "F")

        blk_count += 1

for sIdy in range(nsy):
    for sIdx in range(nsx):

        index = sIdx + sIdy*nsx

        SPS[index, 0] = sx[sIdy, sIdx]
        SPS[index, 1] = sy[sIdy, sIdx] 
        SPS[index, 2] = 12.5

        actives_x = slice(sIdx, sIdx + spread_x) 
        actives_y = slice(sIdy*spread_y, sIdy*spread_y + spread_y)

        XPS[index, 0] = index
        XPS[index, 1] = int(np.min(rec_index[actives_y, actives_x]))
        XPS[index, 2] = int(np.max(rec_index[actives_y, actives_x])) + 1

np.savetxt("../inputs/geometry/migration_test_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/migration_test_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/migration_test_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")
