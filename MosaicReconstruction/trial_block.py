# Reconstruction of stitched projections

import tomopy
import time
import h5py
import numpy as np
import os

#### parameter to vary ####
slice_range = [1024, 1024 + 32]

method='gridrec'     # 'fbp' takes very long
iterK = 2            # number of iterations for iterative methods
pad_sinogram = 3000

ncore = 4

###### Input: Dataset to process
paramfile = '/home/mattia/Documents/Cerebellum22/MosaicReconstruction/example/param_files/cerebellum_tile6.txt'
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
pixsize_um = float(infoStruct['pixsize_um'])
outputcrop = infoStruct['outputcrop'].split(",")
for i in range(len(outputcrop)):
    outputcrop[i] = int(outputcrop[i])

###### 0.4 Set up directories
projdir = projbasedir + samplename + os.sep + 'proj' + os.sep;
recodir = recobasedir + samplename + os.sep + 'reco' + os.sep;

os.makedirs(recodir, mode=0o755, exist_ok=True)

print('Writing results to ' + recodir)    

###### 0.5 Load measurement info
f = h5py.File(projdir + 'angles.h5','r')
angles = np.pi*np.squeeze(np.array(f['angles']))/180.0
ip180 = angles.shape[0]

pixsize = pixsize_um*1e-6;      # [m]
pixsize_mm = pixsize_um*1e-3;

# info on projections
f = h5py.File(projdir + 'proj_f_' + '%04d' % (1) + '.h5', 'r')   # read information
sx = f['/proj'].shape[0]
sy = f['/proj'].shape[1]
sz = slice_range[1] - slice_range[0]

typ = f['/proj'].dtype

print('Loading slices ' + str(slice_range[0] + 1) + ' to '
	+ str(slice_range[1]) + '...')
t1 = time.time()

# load projection block 
projblock = np.empty((ip180, sz, sx), typ) # theta z x (projection order)
print('Loading projections...')
t2 = time.time()
for p in range(ip180):
	f = h5py.File(projdir + 'proj_f_%04d' % (p+1) +  '.h5', 'r')
	typ = f['/proj'].dtype
	projblock[p,:,:] = f['/proj'][:,slice_range[0]:slice_range[1]]\
		.transpose()

executionTime = (time.time() - t2)
print('read projections: %.2f sec ' % (executionTime))

print('Reconstruction and writing...')
t3 = time.time()
projblock = tomopy.minus_log(projblock)
if pad_sinogram != 0:
	projblock = np.pad(projblock,
		((0,0),(0,0),(pad_sinogram,pad_sinogram)),
		mode='edge')

# sinogram_order means: stack of sinograms (z theta x)
# non-sinogram_order means: stack of projections (theta z x)
if method=='fbp' or method=='gridrec':
	reco = tomopy.recon(projblock, angles,
		sinogram_order=False,algorithm=method,ncore=ncore)\
		/pixsize_mm
else:
	reco = tomopy.recon(projblock, angles,
		sinogram_order=False,algorithm=method,
		num_iter=iterK,ncore=ncore)/pixsize_mm
    
# remove the padding
reco = reco[:,pad_sinogram:-pad_sinogram,pad_sinogram:-pad_sinogram]
# crop
reco = reco[:,outputcrop[2]:sx-outputcrop[3],outputcrop[0]:sx-outputcrop[1]];

# write HDF5
outfname = '%sreco_%05d-%05d.h5' % ((recodir,slice_range[0] + 1,
	slice_range[1]))
fid = h5py.File(outfname, 'w')
ds = fid.create_dataset('/reco', reco.shape, dtype=reco.dtype)
ds[()] = reco
fid.close()

executionTime = (time.time() - t3)
print('reconstruction: %.2f sec ' % (executionTime))

executionTime = (time.time() - t1)
print('total time: %.2f sec ' % (executionTime))
