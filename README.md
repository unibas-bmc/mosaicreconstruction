# Software Pipeline for Extended Field-of-View X-ray Tomography by Projection Stitching


## Introduction

To extend the imaged field-of-view, projections are acquired with
offset rotation axis and vertical stage position. The projections
are stitched laterally and vertically before tomographic
reconstruction.

This software pipeline implements the stitching of projections
in the lateral (x) and vertical (y) direction.
The stitching positions are estimated automatically by
maximizing image cross-correlation, starting from the motor
positions. The user has the option to manually verify and adjust the
results based on trial reconstructions. We use
[TomoPy](https://tomopy.readthedocs.io) for tomographic reconstruction.

The organisation of the raw data is defined by an *info file*
(xlsx or csv format), and the pipeline is configured by the
*param file* (plain text format).


## Info File

The info file (e.g. [example/param_files/info.csv](
example/param_files/info.csv)) lists all scan directories, and to which
heightstep and ring the data belong.
The rings are numbered 1 to *n*, starting from the inside and the
heightsteps 1 to *m*, starting form the top.
The scan names are constructed from the `scandir` and the `suffix`,
e.g. `001_sample123_ring` and `z1_x2` would be combined
to `001_sample123_ring_z1_x2`.


## Param File

The param file (e.g. [example/param_files/param.txt](
example/param_files/param.txt)) contains configuration
settings for stitching and reconstruction. These include
 * sample name
 * path to the info file
 * paths to raw and processed data
 * parameters for Paganin or Gaussian filtering (optional)
 * scaling of gray values in the reconstructed slices


# Main Steps

1. Determine overlap positions for lateral stitching:
	```
	OverlapFinderX(paramfile);
	```

2. Perform lateral stitching (stitched projections are
	saved in `[projpath "stitched_proj_filtered"]`):
	```
	ProjectionProcessing_pass1(paramfile, manstitchposx);
	```

3. Determine overlap positions for vertical stitching:
	```
	OverlapFinderY(paramfile);
	```

4. Perform vertical stitching (stitched projections are
	saved in `[projpath2 "proj"]`):
	```
	ProjectionProcessing_pass2(paramfile, manstitchposy);
	```

5. Tomographic reconstruction (reconstructed slices are
	saved in `[recopath "reco"]`):
	```
	python utils/block_reconstruction.py example/param_files/param.txt
	```


## Manual intervention

Especially for large specimens without highly-attenuating structures,
contrast-to-noise ratio is typically very low in the projections.
This makes the automated search for overlap positions more
error-prone, and the need for manual intervention more likely.
The Matlab script [MosaicReco.m](MosaicReco.m) contains code snippets
to visually verify and manually optimize stitching and reconstruction
parameters. It is not intended to be run as a whole, but instead
divided into sections to be run individually.


## Dependencies

Requires the python packages
[TomoPy](https://tomopy.readthedocs.io),
[h5py](https://www.h5py.org/), and
[Tifffile](https://github.com/cgohlke/tifffile/).
