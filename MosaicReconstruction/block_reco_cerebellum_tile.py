# Reconstruction of stitched projections

#import code       # for debugging, breakpoint with
                   # code.interact(local=locals())
import time
import os
import tomopy
import h5py
import numpy as np
import pandas

import libtiff as tif
import threading
import gc


###### settings, should all come from parameter file
doBlock = False         # True: to block of 32 slices at once,
                        # else slice-wise reconstruction
method='gridrec'        # 'fbp' takes very long
iterK = 2               # number of iterations for iterative methods
verboseExtra = True     # True: more information
pad_sinogram = 3000	# pad sinogram on both sides with n pixels

###### Input: Dataset to process
paramfile = "/home/mattia/Documents/Cerebellum22/MosaicReconstruction/example/param_files/cerebellum_tile3.txt"
print('Using ' + paramfile)

###### 0.2 Read param file
# read file, make structure array
infoStruct = {}
with open(paramfile,'r') as f:
    for line in f: # iterate over each line
        allSplit = line.split("\t")
        if len(allSplit)==2:
            # expecting two values with a tab separated
            name, valueAndMore = line.split("\t") # split it by tabs
            value, others = valueAndMore.split("\n")  # remove endline
            infoStruct[name]=value

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
ncore = int(infoStruct['ncore'])
pixsize_um = float(infoStruct['pixsize_um'])

###### 0.4 Set up directories
projdir = projbasedir + samplename + os.sep + 'proj' + os.sep
recodir = recobasedir + samplename + os.sep + 'reco' + os.sep

print('Writing results to ' + recodir)

###### 0.5 Load measurement info
f = h5py.File(projdir + 'angles.h5','r')
angles = np.pi*np.squeeze(np.array(f['angles']))/180.0
ip180 = angles.shape[0]

T = pandas.read_csv(infofile)
pixsize = pixsize_um*1e-6       # [m]
pixsize_mm = pixsize_um*1e-3

# info on projections
f = h5py.File(projdir + 'proj_f_' + '%04d' % (1) + '.h5', 'r')
sx = f['/proj'].shape[0]
sy = f['/proj'].shape[1]
nblocks = int(sy/blocksize)

width = sx
height = sy
typ = f['/proj'].dtype
f.close()
del f

def load_projection_block(projdir, b, blocksize, sx, ip180, typ):
    tyi = b*blocksize+1    # this y initial (indices 1 based)
    tyf = (b+1)*blocksize  # this y final
    projblock = np.empty((blocksize, sx, ip180), typ)
    print('Loading projections...')
    t2 = time.time()
    for p in range(ip180):
        f = h5py.File(projdir + 'proj_f_%04d' % (p+1) +  '.h5', 'r')
        projblock[:,:,p] = f['/proj'][:,tyi-1:tyf].T
        f.close()
        del f
    executionTime = (time.time() - t2)
    print('read projections: %.2f sec ' % (executionTime))
    return projblock

class LoadBlockThread(threading.Thread):
    def __init__(self, *args):
        self.args = args
        self.projblock = None
        super(LoadBlockThread, self).__init__()
    def run(self):
        self.projblock = load_projection_block(*self.args)

print('Loading block 1/' + str(nblocks) + '...')

load_block_thread = LoadBlockThread(projdir, 0, blocksize,
    sx, ip180, typ)
load_block_thread.start()

###### 1.1 loop over y, load singogram, run reco, output reconstruction
for b in range(nblocks):

    tyi = b*blocksize+1    # this y initial (indices 1 based)
    tyf = (b+1)*blocksize  # this y final
    tyr = range(tyi,tyf+1)   # this y range

    # retrieve projection block
    load_block_thread.join()
    projblock = load_block_thread.projblock

    # load next projection block
    if b + 1 < nblocks:
        load_block_thread = LoadBlockThread(projdir, b + 1, blocksize,
            sx, ip180, typ)
        load_block_thread.start()

    print('Reconstructing block ' + str(b+1) + '/' + str(nblocks) + '...')
    # reshape into sinogram stack
    projblock = np.transpose(projblock,(1,2,0))
    sz=projblock.shape

    print('Reconstruction and writing...')
    t3 = time.time()
    if doBlock:
        # do strip at once
        this_sino_log = tomopy.minus_log(projblock, ncore=ncore)
        this_sino_log = np.transpose(this_sino_log,(2,1,0))
        if pad_sinogram != 0:
            this_sino_log = np.pad(this_sino_log,
                ((0,0),(0,0),(pad_sinogram,pad_sinogram)),
                mode='edge')

        if method=='fbp' or method=='gridrec':
            reco = tomopy.recon(this_sino_log, angles,
                sinogram_order=True, algorithm=method,
                ncore=ncore)/pixsize_mm
        else:
            reco = tomopy.recon(this_sino_log, angles,
                sinogram_order=True, algorithm=method,
                num_iter=iterK, ncore=ncore)/pixsize_mm

        # remove the padding
        if pad_sinogram != 0:
            reco = reco[:,pad_sinogram:-pad_sinogram,
                pad_sinogram:-pad_sinogram]
        # crop
        reco = reco[:,outputcrop[2]:sz[0]-outputcrop[3],
            outputcrop[0]:sz[0]-outputcrop[1]]
        # convert to uint16
        # reco = np.uint16(2**16*((reco-outputgrayscale[0])
        #     /(outputgrayscale[1]-outputgrayscale[0])))
        reco = np.int16((2**16*((reco-outputgrayscale[0])
            /(outputgrayscale[1]-outputgrayscale[0])))-2**15)

        for iy in range(len(tyr)):
            outfname = '%s/reco_%05d.tif' % ((recodir,tyr[iy]))
            fid = tif.TIFF.open(outfname, 'w')
            fid.write_image(np.squeeze(reco[iy,:,:]))
            fid.close()
        del reco, this_sino_log
    else:
        for iy in range(len(tyr)):
            this_sino_log = tomopy.minus_log(projblock[:,:,iy])
            this_sino_log = np.transpose(this_sino_log)
            this_sino_log = np.expand_dims(this_sino_log,axis=0)
            if pad_sinogram != 0:
                this_sino_log = np.pad(this_sino_log,
                    ((0,0),(0,0),(pad_sinogram,pad_sinogram)),
                    mode='edge')

            if method=='fbp' or method=='gridrec':
                reco = tomopy.recon(this_sino_log, angles,
                sinogram_order=True,algorithm=method)/pixsize_mm
            else:
                reco = tomopy.recon(this_sino_log, angles,
                    sinogram_order=True, algorithm=method,
                    num_iter=iterK)/pixsize_mm

            # remove the padding
            if pad_sinogram != 0:
                reco = reco[:,pad_sinogram:-pad_sinogram,
                    pad_sinogram:-pad_sinogram]
            reco = np.squeeze(reco)
            # crop
            reco = reco[outputcrop[2]:sz[0]-outputcrop[3],
                outputcrop[0]:sz[0]-outputcrop[1]]
            # convert to int16
            # reco = np.uint16(2**16*((reco-outputgrayscale[0])
            #     /(outputgrayscale[1]-outputgrayscale[0])))
            reco = np.int16((2**16*((reco-outputgrayscale[0])
                /(outputgrayscale[1]-outputgrayscale[0])))-2**15)

            outfname = '%s/reco_%05d.tif' % ((recodir,tyr[iy]))
            fid = tif.TIFF.open(outfname, 'w')
            fid.write_image(reco)
            fid.close()
            del reco, this_sino_log

    del projblock
    executionTime = (time.time() - t3)
    print('reconstruction: %.2f sec ' % (executionTime))
    gc.collect()
