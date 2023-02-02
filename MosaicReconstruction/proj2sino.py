# Reconstruction of stitched projections

import time
import os
import tomopy
import h5py
import numpy as np

import threading
import gc

t0 = time.time()

# Input: Dataset to process
paramfile = "/home/mattia/Documents/Cerebellum22/MosaicReconstruction/example/param_files/cerebellum_tile7.txt"
print('Using ' + paramfile)

# 0.2 Read param file
# read file, make structure array
infoStruct = {}
with open(paramfile, 'r') as f:
    for line in f:  # iterate over each line
        allSplit = line.split("\t")
        if len(allSplit) == 2:
            # expecting two values with a tab separated
            name, valueAndMore = line.split("\t")  # split it by tabs
            value, others = valueAndMore.split("\n")  # remove endline
            infoStruct[name] = value

f.close()
# assign variables
samplename = infoStruct['samplename']
projbasedir = infoStruct['projpath2']
recobasedir = infoStruct['recopath']
stripheight = float(infoStruct['stripheight'])

## Overwrite blocksize to use memory more efficiently, as we're
## not computing anything
# blocksize = int(infoStruct['stripheight'])
blocksize = 64

# 0.4 Set up directories
projdir = projbasedir + samplename + os.sep + 'proj' + os.sep
sinodir = projbasedir + samplename + os.sep + 'sino' + os.sep

os.makedirs(sinodir, mode=0o755, exist_ok=True)

print('Writing results to ' + sinodir)

# 0.5 Load measurement info
f = h5py.File(projdir + 'angles.h5', 'r')
ip180 = f['angles'].shape[1]
f.close()

# info on projections
f = h5py.File(projdir + 'proj_f_' + '%04d' % 1 + '.h5', 'r')
sx = f['/proj'].shape[0]
sy = f['/proj'].shape[1]
nblocks = int(sy / blocksize)

typ = f['/proj'].dtype
f.close()
del f


def load_projection_block(_projdir, _b, _blocksize, _sx, _ip180, _typ):
    _tyi = _b * _blocksize + 1  # this y initial (indices 1 based)
    _tyf = (_b + 1) * _blocksize  # this y final
    _projblock = np.empty((_ip180, _blocksize, _sx), _typ)
    for p in range(_ip180):
        _f = h5py.File(_projdir + 'proj_f_%04d' % (p + 1) + '.h5', 'r')
        _projblock[p, :, :] = _f['/proj'][:, _tyi - 1:_tyf].T
        _f.close()
        del _f
    return _projblock


# 1.1 loop over y, build and store singogram
for b in range(nblocks):

    print("Processing block {}/{}".format(b+1, nblocks))

    tyi = b * blocksize + 1  # this y initial (indices 1 based)
    tyf = (b + 1) * blocksize  # this y final
    tyr = range(tyi, tyf + 1)  # this y range

    t1 = time.time()
    # retrieve projection block
    projblock = load_projection_block(
            projdir, b, blocksize, sx, ip180, typ)
    assert(projblock.shape == (ip180, blocksize, sx))

    t2 = time.time()
    print("Loading time: {} s".format(t2 - t1))

    # reshape to synogram stack
    projblock = np.transpose(projblock, (1, 0, 2))
    assert(projblock.shape == (blocksize, ip180, sx))

    for iy in range(len(tyr)):
        outfname = '%ssino_%05d.h5' % (sinodir, tyr[iy])
        f = h5py.File(outfname, 'w')
        f.create_dataset("/sino", data=projblock[iy,:,:])
        f.close()
    del projblock
    t3 = time.time()
    print("Storing time: {} s".format(t3 - t2))

    gc.collect()

execution_time = time.time() - t0
print("Total time: {} s".format(execution_time))
