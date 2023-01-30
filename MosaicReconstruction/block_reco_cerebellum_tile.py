# Reconstruction of stitched projections

import time
import os
import tomopy
import h5py
import numpy as np

import libtiff as tif
import threading
import gc

# settings, should all come from parameter file

ncore = 4

method = 'gridrec'  # 'fbp' takes very long
iterK = 2  # number of iterations for iterative methods
verboseExtra = True  # True: more information
pad_sinogram = 3000  # pad sinogram on both sides with n pixels

# Input: Dataset to process
paramfile = "/home/mattia/Documents/Cerebellum22/MosaicReconstruction/example/param_files/cerebellum_tile5.txt"
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
infofile = infoStruct['infofile']
projbasedir = infoStruct['projpath2']
recobasedir = infoStruct['recopath']
stripheight = float(infoStruct['stripheight'])
outputcrop = infoStruct['outputcrop'].split(",")
for i in range(len(outputcrop)):
    outputcrop[i] = int(outputcrop[i])

outputgrayscale = infoStruct['outputgrayscale'].split(",")
for i in range(len(outputgrayscale)):
    outputgrayscale[i] = float(outputgrayscale[i])

blocksize = int(infoStruct['stripheight'])
# ncore = int(infoStruct['ncore'])
pixsize_um = float(infoStruct['pixsize_um'])

# 0.4 Set up directories
projdir = projbasedir + samplename + os.sep + 'proj' + os.sep
sinodir = projbasedir + samplename + os.sep + 'sino' + os.sep
recodir = recobasedir + samplename + os.sep + 'reco' + os.sep

print('Writing results to ' + recodir)

# 0.5 Load measurement info
f = h5py.File(projdir + 'angles.h5', 'r')
angles = np.pi * np.squeeze(np.array(f['angles'])) / 180.0
ip180 = angles.shape[0]

pixsize = pixsize_um * 1e-6  # [m]
pixsize_mm = pixsize_um * 1e-3


# info on sinograms
f = h5py.File(sinodir + 'sino_' + '%05d' % (1) + '.h5', 'r')
assert(f['/sino'].shape[0] == ip180)
sx = f['/sino'].shape[1]
typ = f['/sino'].dtype
f.close()

# info on projections
f = h5py.File(projdir + 'proj_f_' + '%04d' % (1) + '.h5', 'r')
assert(f['/proj'].shape[0] == sx)
sz = f['/proj'].shape[1]

nblocks = int(sz / ncore)

f.close()
del f

# 1.1 loop over y, load singogram, run reco, output reconstruction
for b in range(nblocks):

    tyi = b * ncore + 1  # this y initial (indices 1 based)
    tyf = (b + 1) * ncore  # this y final
    tyr = range(tyi, tyf + 1)  # this y range

    sinoblock = np.empty((ncore, ip180, sx), typ)
    print('Loading sinograms...')
    t2 = time.time()
    for s in range(ncore):
        f = h5py.File(sinodir + 'sino_%05d' % (tyr[s]) + '.h5', 'r')
        sinoblock[s, :, :] = f['/sino'][()]
        f.close()
        del f
    execution_time = (time.time() - t2)
    print('read sinograms: %.2f sec ' % execution_time)

    print('Reconstructing block ' + str(b + 1) + '/' + str(nblocks) + '...')

    print('Reconstruction and writing...')
    t3 = time.time()

    sinoblock = tomopy.minus_log(sinoblock, ncore=ncore)
    sinoblock = np.transpose(sinoblock, (2, 1, 0))
    if pad_sinogram != 0:
        sinoblock = np.pad(sinoblock,
                               ((0, 0), (0, 0), (pad_sinogram, pad_sinogram)),
                               mode='edge')

    if method == 'fbp' or method == 'gridrec':
        reco = tomopy.recon(sinoblock, angles,
                            sinogram_order=True, algorithm=method,
                            ncore=ncore) / pixsize_mm
    else:
        reco = tomopy.recon(sinoblock, angles,
                            sinogram_order=True, algorithm=method,
                            num_iter=iterK, ncore=ncore) / pixsize_mm

    # remove the padding
    if pad_sinogram != 0:
        reco = reco[:, pad_sinogram:-pad_sinogram,
                pad_sinogram:-pad_sinogram]
    # crop
    reco = reco[:, outputcrop[2]:-outputcrop[3],
            outputcrop[0]:-outputcrop[1]]
    # convert to uint16
    # reco = np.uint16(2**16*((reco-outputgrayscale[0])
    #     /(outputgrayscale[1]-outputgrayscale[0])))
    reco = np.int16((2 ** 16 * ((reco - outputgrayscale[0])
                                / (outputgrayscale[1] - outputgrayscale[0]))) - 2 ** 15)

    for iy in range(len(tyr)):
        outfname = '%s/reco_%05d.tif' % (recodir, tyr[iy])
        fid = tif.TIFF.open(outfname, 'w')
        fid.write_image(np.squeeze(reco[iy, :, :]))
        fid.close()
    del reco, sinoblock

    executionTime = (time.time() - t3)
    print('reconstruction: %.2f sec ' % executionTime)
    gc.collect()
