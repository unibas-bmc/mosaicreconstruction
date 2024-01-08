#!/usr/bin/env python
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt

rawdir = "/media/mattia/Soleil_2023_Disk2/nxs"

orig = "032_Mouse89883_5823_eth_9000proj_ring_z6_x1"
# orig = "032_Mouse89883_5823_eth_9000proj_ring_redo_z5_x3"

# dummy = "032_Mouse89883_5823_eth_9000proj_ring_z6_x2_dummy"
# r0 = -287.0 / 2. + 2 * 1759.4
# r1 = -287.0 / 2. + 2 * 1759.4 + 2048
dummy = "032_Mouse89883_5823_eth_9000proj_ring_z6_x3_dummy"
r0 = -287.0 / 2. + 3 * 1759.4
r1 = -287.0 / 2. + 3 * 1759.4 + 2048

#### CYLINDER PROFILE ####
ra = -287.0 / 2. + 1759.4 + 2048.
rb = 7179.
A = 0.79957345
B = 0.8962663

C = np.log(A) / np.log(B)

R = np.sqrt((C ** 2 * rb ** 2 - ra ** 2) / (C ** 2 - 1.))
mu = np.abs(np.log(B) / (2 * np.sqrt(R ** 2 - rb ** 2)))
assert(np.allclose(mu, np.abs(np.log(A)
    / (2 * np.sqrt(R ** 2 - ra ** 2)))))

print("R =", R)
print("mu =", mu)

def T(x, R, mu):
    return np.exp(-2. * mu * R * np.sqrt(1. - (x / R) ** 2))

assert(np.allclose(T(ra, R, mu), A))
assert(np.allclose(T(rb, R, mu), B))
#### /CYLINDER PROFILE ####

# load and flat-correct orig
if False:
    ref = np.mean(h5py.File(os.path.join(rawdir, orig, "pre_ref.nxs"))\
            ["/ref_0"][:,:,-1].astype(np.float32), axis=0)
    ref += np.mean(h5py.File(os.path.join(rawdir, orig, "post_ref.nxs"))\
            ["/ref_0"][:,:,-1].astype(np.float32), axis=0)
    ref /= 2.
    dark = np.mean(h5py.File(os.path.join(rawdir, orig, "post_dark.nxs"))\
            ["/dark_0"][:,:,-1].astype(np.float32), axis=0)
    im = h5py.File(os.path.join(rawdir, orig, orig + ".nxs"))\
            ["/flyscan_0003/scan_data/orca_image"][:,:,-1]\
            .astype(np.float32)
    im = (im - dark) / (ref - dark)
    plt.imshow(im)
    plt.show()

# alter dummy
if True:
    r = np.linspace(r0, r1 - 1., 2048)
    s = T(r, R, mu).reshape((1, 1, -1))
    u16_max = 65535
    fill_value = np.round(s * u16_max).astype(np.uint16)
    f = h5py.File(os.path.join(rawdir, dummy, "pre_ref.nxs"), "r+")
    f["/ref_0"][()] = u16_max
    f.close()
    f = h5py.File(os.path.join(rawdir, dummy, "post_ref.nxs"), "r+")
    f["/ref_0"][()] = u16_max
    f.close()
    f = h5py.File(os.path.join(rawdir, dummy, "post_dark.nxs"), "r+")
    f["/dark_0"][()] = 0
    f.close()
    f = h5py.File(os.path.join(rawdir, dummy, dummy + ".nxs"), "r+")
    f["/flyscan_0003/scan_data/orca_image"][()] = fill_value
    f.close()
