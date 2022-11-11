# Reconstruction of stitched projections

#import code       # for debugging, breakpoint with code.interact(local=locals())
import tomopy
import pylab
import time
import h5py
import numpy as np
import os
import pandas
import sys

import libtiff as tif

#### parameter to vary ####
slice_no = 128
outputgrayscale = [0.0, 0.12]

method='gridrec'     # 'fbp' takes very long
iterK = 2;           # number of iterations for iterative methods
pad_sinogram = 3000

###### Input: Dataset to process
paramfile = './example/param_files/Cerebellum_8Rings.txt'
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

###### 0.4 Set up directories
projdir = projbasedir + samplename + os.sep + 'proj' + os.sep;
recodir = recobasedir + samplename + os.sep + 'reco' + os.sep;

os.makedirs(recodir, mode=0o755, exist_ok=True)

print('Writing results to ' + recodir)    

###### 0.5 Load measurement info
f = h5py.File(projdir + 'angles.h5','r')
angles = np.pi*np.squeeze(np.array(f['angles']))/180.0
ip180 = angles.shape[0]

T = pandas.read_excel(infofile)
pixsize_um = 0.65 # TODO: read from par file
pixsize = pixsize_um*1e-6;      # [m]
pixsize_mm = pixsize_um*1e-3;

# info on projections
t = tif.TIFF.open(projdir + 'proj_f_' + '%04d' % (1) + '.tif')   # read information
sx = t.GetField("ImageWidth")
sy = t.GetField("ImageLength")
blocksize = t.GetField("RowsPerStrip")
nblocks = int(sy/blocksize)

width = t.GetField("ImageWidth")
height = t.GetField("ImageLength")
bits = t.GetField('BitsPerSample')
sample_format = t.GetField('SampleFormat')
typ = t.get_numpy_type(bits, sample_format)

b = int(slice_no / blocksize)
print('Loading block ' + str(b+1) + '/' + str(nblocks) + '...')
t1 = time.time()
tyi = b*blocksize+1;   # this y initial
tyf = (b+1)*blocksize; # this y final
tyr = range(tyi,tyf+1);  # this y range

# load projection block 
projblock = np.empty((blocksize, sx, ip180), typ)
print('Loading projections...')
t2 = time.time()
for p in range(ip180):
	t = tif.TIFF.open(projdir + 'proj_f_%04d' % (p+1) +  '.tif', 'r')
	bits = t.GetField('BitsPerSample')
	sample_format = t.GetField('SampleFormat')
	typ = t.get_numpy_type(bits, sample_format)
	tmp = np.empty((blocksize, width), typ)
	size = tmp.nbytes
	ReadStrip = t.ReadEncodedStrip
	elem = ReadStrip(b, tmp.ctypes.data, size)
	#tmp[np.isnan(tmp)] = 1;
	projblock[:,:,p] = tmp;

executionTime = (time.time() - t2)
print('read projections: %.2f sec ' % (executionTime))

# reshape into sinogram stack
projblock = np.transpose(projblock,(1,2,0));
sz=projblock.shape

print('Reconstruction and writing...')
t3 = time.time()
iy = slice_no % blocksize
this_sino_log = tomopy.minus_log(projblock[:,:,iy])
del projblock
this_sino_log = np.transpose(this_sino_log)
this_sino_log = np.expand_dims(this_sino_log,axis=0)
if pad_sinogram != 0:
	this_sino_log = np.pad(this_sino_log,
		((0,0),(0,0),(pad_sinogram,pad_sinogram)),
		mode='edge')

if method=='fbp' or method=='gridrec':
	reco = tomopy.recon(this_sino_log, angles,
		sinogram_order='False',algorithm=method)/pixsize_mm
else:
	reco = tomopy.recon(this_sino_log, angles,
		sinogram_order='False',algorithm=method,
		num_iter=iterK)/pixsize_mm
    
reco = np.squeeze(reco)
# remove the padding
reco = reco[pad_sinogram:-pad_sinogram,pad_sinogram:-pad_sinogram]
# crop
reco = reco[outputcrop[2]:sz[0]-outputcrop[3],outputcrop[0]:sz[0]-outputcrop[1]];

# write HDF5
outfname = '%sreco_%05d.h5' % ((recodir,tyr[iy]))
fid = h5py.File(outfname, 'w')
ds = fid.create_dataset('/reco', reco.shape, dtype=reco.dtype)
ds[()] = reco
fid.close()

# convert to int16
# reco = np.uint16(2**16*((reco-outputgrayscale[0])/(outputgrayscale[1]-outputgrayscale[0])));
reco = np.int16((2**16*((reco-outputgrayscale[0])/(outputgrayscale[1]-outputgrayscale[0])))-2**15);

outfname = '%sreco_%05d.tif' % ((recodir,tyr[iy]))
fid = tif.TIFF.open(outfname, 'w')
fid.write_image(reco)
fid.close()

executionTime = (time.time() - t3)
print('reconstruction: %.2f sec ' % (executionTime))

executionTime = (time.time() - t1)
print('total time: %.2f sec ' % (executionTime))
