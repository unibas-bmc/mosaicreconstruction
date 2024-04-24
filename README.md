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
