#!/usr/bin/env python

import tomopy
import numpy as np
import SimpleITK as sitk
import sys


ratio = 0.6

fixed_file = "/media/griffin/griffin_disk11/anatomix_mousebrain/" \
    + "mouse6_perf_eth/stitchpos_tests/overlap1_reco.mha"
filename = "/media/griffin/griffin_disk11/anatomix_mousebrain/" \
    + "mouse6_perf_eth/stitchpos_tests/overlap_mask.mha"

if len(sys.argv) > 1:
    fixed_file = sys.argv[1]
if len(sys.argv) > 2:
    filename = sys.argv[2]
if len(sys.argv) > 3:
    ratio = float(sys.argv[3])

fixed = sitk.ReadImage(fixed_file)
sx, sy, sz = fixed.GetSize()

mask = np.ones((sz, sy, sx))
mask = tomopy.circ_mask(mask, axis=0, ratio=ratio).astype(np.uint8)

mask_sitk = sitk.GetImageFromArray(mask)
mask_sitk.CopyInformation(fixed)
sitk.WriteImage(mask_sitk, filename)

print(filename)
