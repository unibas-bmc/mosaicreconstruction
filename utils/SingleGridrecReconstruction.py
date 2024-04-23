#import code       # for debugging, breakpoint with code.interact(local=locals())
import tomopy
import time
import h5py
import numpy as np
#import libtiff as tif
import os
import sys

#import matplotlib.pyplot as plt
#%matplotlib qt
#import pandas
#import pylab

###### settings, should all come from parameter file
method='gridrec'     # 'fbp' takes very long
iterK = 2;           # number of iterations for iterative methods

###### inputs: sinogram file
#sinogramfile = '/media/griffin/griffin_disk7/soleil_feb2021/20200712_mueller_processed/brain4_ethanol/stitching_testrecos/test_sino.h5'
sinogramfile = sys.argv[1]
print('Using ' + sinogramfile)

# read sizes, angles, and .par file info
hf = h5py.File(sinogramfile, 'a')
n1 = hf.get('/angles')
angles = np.array(n1)
n2 = hf.get('/sino')
sino = np.array(n2)
n3 = hf.get('/pixsize_mm')
pixsize_mm = np.array(n3)
hf.close()

# convert angle degrees to radians
angles = np.pi*angles/180.0

###### run reco, output reconstruction
t3 = time.time()
this_sino_log = tomopy.minus_log(sino)
if sino.ndim<3:
    this_sino_log = np.expand_dims(this_sino_log,axis=0)
    
if method=='fbp' or method=='gridrec':
    reco = tomopy.recon(this_sino_log, angles, sinogram_order='False',algorithm=method)/pixsize_mm
else:
    reco = tomopy.recon(this_sino_log, angles, sinogram_order='False',algorithm=method,num_iter=iterK)/pixsize_mm
reco = np.squeeze(reco,axis=0)

#addstr = '_reco.tif'
#outfullname = os.path.splitext(sinogramfile)[0] + addstr
#fid = tif.TIFF.open(outfullname, 'w')
#fid.write_image(reco)
#fid.close()

addstr = '_reco.h5'
outfullname = os.path.splitext(sinogramfile)[0] + addstr
hf = h5py.File(outfullname, 'w')
hf.create_dataset('reco', data=reco)
hf.close()

executionTime = (time.time() - t3)
print('reconstruction time: %.2f sec ' % (executionTime))
    
        
