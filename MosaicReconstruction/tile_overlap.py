#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import norm
import pandas
import tifffile
import SimpleITK as sitk
import os


class RecoParam:
    def __init__(self, paramfile):
        print("Reading parameter file: {}".format(paramfile))
        # read file, make structure array
        infoStruct = {}
        with open(paramfile, 'r') as f:
            for line in f:  # iterate over each line
                allSplit = line.split("\t")
                if len(allSplit) == 2:
                    # expecting two values with a tab separated
                    # split it by tabs
                    name, valueAndMore = line.split("\t")
                    # remove endline
                    value, others = valueAndMore.split("\n")
                    infoStruct[name] = value
            f.close()

        # assign variables
        self.samplename = infoStruct['samplename']
        self.infofile = infoStruct['infofile']
        self.rawdatapath = infoStruct['rawdatapath']
        self.projbasedir = infoStruct['projpath2']
        self.recobasedir = infoStruct['recopath']
        self.stripheight = float(infoStruct['stripheight'])
        self.outputcrop = infoStruct['outputcrop'].split(",")
        for i in range(len(self.outputcrop)):
            self.outputcrop[i] = int(self.outputcrop[i])

        self.outputgrayscale = infoStruct['outputgrayscale'].split(",")
        for i in range(len(self.outputgrayscale)):
            self.outputgrayscale[i] = float(self.outputgrayscale[i])

        self.blocksize = int(infoStruct['stripheight'])
        self.ncore = int(infoStruct['ncore'])
        self.pixsize_um = float(infoStruct['pixsize_um'])
        print("Reading info file: {}".format(self.infofile))
        self.infotable = pandas.read_csv(self.infofile)

        ## read motor positions
        scan_center = self.infotable[pandas.DataFrame([
                self.infotable.ring == 1,
                self.infotable.heightstep == 1,
            ]).all()]
        assert(scan_center.shape[0] == 1)
        scan_center = scan_center.iloc[0]
        scanname = scan_center.scandir
        if scan_center.suffix != 0:
            scanname += "_{}".format(scan_center.suffix)
        parfile = "{}{}/{}.par".format(self.rawdatapath,
            scanname, scanname)

        print("Reading PyHST .par file: {}".format(parfile))
        with open(parfile, 'r') as f:
            for l in f:
                l = l.split(' ')
                if l[0] == "#tomo1tx":
                    assert(l[1] == "=")
                    self.tomo1tx = float(l[2])
                elif l[0] == "#tomo1tz":
                    assert(l[1] == "=")
                    self.tomo1tz = float(l[2])
                elif l[0] == "#tomo1tu1":
                    assert(l[1] == "=")
                    self.tomo1tu1 = float(l[2])
                elif l[0] == "#tomo1tv1":
                    assert(l[1] == "=")
                    self.tomo1tv1 = float(l[2])

        # read reco width
        recopath = self.recobasedir + self.samplename + "/reco/"
        print("Reading slice from: {}".format(recopath))
        tif = tifffile.TiffFile(recopath + "reco_00001.tif")
        assert(tif.pages[0].shape[0] == tif.pages[0].shape[1])
        self.recowidth = tif.pages[0].shape[0]


def rect_is_contained(r, x, y, R):
    d = norm(r - x.reshape(1,2), axis=1)
    if not np.all(d < R):
        return False
    d = norm(r - y.reshape(1,2), axis=1)
    if not np.all(d < R):
        return False
    return True


def rect_grow(r):
    return r + np.array([
            [-1, -1],
            [+1, -1],
            [+1, +1],
            [-1, +1],
        ])


def overlap_rects(param1, param2, ratio=0.9, plot=False):

    print("Find overlap for radii with ratio: {}".format(ratio))
    width = param1.recowidth # pixels
    radius = np.round((width + param1.outputcrop[0]
            + param1.outputcrop[1]) / 2)
    print("width: {}".format(width))
    print("radius: {}".format(radius))
    # reduce the radius by a certain factor, so we are out of the range
    # of edge effects and hopefully, to a degree, cupping
    radius *= ratio

    pixsize = param1.pixsize_um * 1e-3

    # u & v motor positions (u is x, v is y)
    x1 = np.array([param1.tomo1tu1, param1.tomo1tv1]) # mm
    x2 = np.array([param2.tomo1tu1, param2.tomo1tv1]) # mm

    x1 /= pixsize
    x2 /= pixsize

    x2 -= x1
    x1 = np.array([np.ceil(width / 2), np.ceil(width/2)])
    x2 += x1
    print("centers {} and {}".format(x1, x2))

    center = np.round((x1 + x2) / 2.)
    rect = np.broadcast_to(center, (4, 2))

    while rect_is_contained(tmp := rect_grow(rect), x1, x2, radius):
        rect = tmp

    assert(rect[1,0] - rect[0,0] >= 200)
    assert(rect[3,1] - rect[0,1] >= 200)

    rect2 = np.round(rect - (x2 - x1).reshape(1,2))
    rect = rect.astype(int)
    print("rect: {}".format(rect))
    rect2 = rect2.astype(int)
    # assume order z, y, x
    # keep in mind that y is stored from top to bottom
    s1 = np.s_[(param1.recowidth-rect[3,1]):(param1.recowidth-rect[0,1]),
            rect[0,0]:rect[1,0]]
    s2 = np.s_[(param2.recowidth-rect2[3,1]):(param2.recowidth-rect2[0,1]),
            rect2[0,0]:rect2[1,0]]

    if plot:
        fig, ax = plt.subplots()
        ax.set_aspect(1)
        ax.set_xlim(10000, 16000)
        ax.set_ylim(0, 8000)

        c1 = plt.Circle(x1, radius, fill=False)
        c2 = plt.Circle(x2, radius, fill=False)
        # c1 = plt.Circle((0.5, 0.5), 0.2, fill=False)
        # c2 = plt.Circle((0.3, 0.7), 0.2, fill=False)
        ax.add_artist(c1)
        ax.add_artist(c2)

        r1 = plt.Rectangle(rect[0,:], rect[1,0] - rect[0,0],
            rect[3,1] - rect[1,1], fill=False)
        ax.add_artist(r1)

        scan1 = plt.Rectangle([0, 0], param1.recowidth, \
                param1.recowidth, fill=False)
        ax.add_artist(scan1)
        scan2 = plt.Rectangle(x2 - x1, param2.recowidth, \
                param2.recowidth, fill=False)
        ax.add_artist(scan2)

        plt.show()

    return s1, s2


def load_volume(param, roi):
    # roi passed as slices in the zyx order
    zrange = [roi[0].start, roi[0].stop]
    yrange = [roi[1].start, roi[1].stop]
    xrange = [roi[2].start, roi[2].stop]
    sz = zrange[1] - zrange[0]
    sy = yrange[1] - yrange[0]
    sx = xrange[1] - xrange[0]
    vol = np.empty((sz, sy, sx))
    print("Loading volume ({}, {}, {}) from: {}{}/reco/".format(
        sz, sy, sx, param.recobasedir, param.samplename))
    for z in range(sz):
        filename = "{}{}/reco/reco_{:05d}.tif".format(
                param.recobasedir, param.samplename, z + zrange[0] + 1)
        im = tifffile.imread(filename)
        vol[z,:,:] = im[yrange[0]:yrange[1], xrange[0]:xrange[1]]
    return vol


def overlap_displacement(param_dir, tile1, tile2, plot=False):

    print("Fixed tile: {}".format(tile1))
    print("Moving tile: {}".format(tile2))
    print("Read parameter files from: {}".format(param_dir))
    param1 = RecoParam("{}cerebellum_tile{}.txt".format(param_dir, tile1))
    param2 = RecoParam("{}cerebellum_tile{}.txt".format(param_dir, tile2))

    basename = "registration_{}-{}".format(tile1, tile2)
    assert(param1.recobasedir == param2.recobasedir)
    write_dir = "{}cerebellum_tile_all/{}/".format(param1.recobasedir,
            basename)
    os.makedirs(write_dir, mode=0o755, exist_ok=True)

    logfile = write_dir + "Log.txt"
    parfile = write_dir + "Parameter.txt"
    transfile = write_dir + "TransformParameter.txt"
    fixedfile = write_dir + "Fixed.mha"
    movingfile = write_dir + "Moving.mha"
    resultfile = write_dir + "Result.mha"

    s1, s2 = overlap_rects(param1, param2, ratio=0.9, plot=plot)
    sz = np.s_[1024-8:1024+8]
    s1 = [sz, s1[1], s1[0]]
    s2 = [sz, s2[1], s2[0]]
    print("Volume in fixed tile: {}".format(s1))
    print("Volume in moving tile: {}".format(s2))
    vol1 = load_volume(param1, s1)
    vol2 = load_volume(param2, s2)
    if plot:
        plt.figure(); plt.imshow(vol1[0,:,:]);
        plt.figure(); plt.imshow(vol2[0,:,:]);
        plt.show()

    print("Write registration results and log to: {}".format(write_dir))
    vol1 = sitk.GetImageFromArray(vol1)
    vol2 = sitk.GetImageFromArray(vol2)
    sitk.WriteImage(vol1, fixedfile)
    sitk.WriteImage(vol2, movingfile)
    parameter_map = sitk.GetDefaultParameterMap("rigid")

    elastix = sitk.ElastixImageFilter()
    elastix.SetFixedImage(vol1)
    elastix.SetMovingImage(vol2)
    elastix.SetParameterMap(parameter_map)
    elastix.WriteParameterFile(parameter_map, parfile)
    # elastix.SetLogFileName(logfile)
    # elastix.SetLogToFile(True)
    # elastix.SetLogToConsole(False)
    elastix.Execute()

    result = elastix.GetResultImage()
    sitk.WriteImage(result, resultfile)

    transform = elastix.GetTransformParameterMap()[0]
    elastix.WriteParameterFile(transform, transfile)
    transform_parameters = np.array(list(transform[
        "TransformParameters"])).astype(np.float32)
    print("TransformParameters: {}".format(transform_parameters))


def _main():

    param_dir = "/home/mattia/Documents/Cerebellum22/" \
            + "MosaicReconstruction/example/param_files/"
    tile1 = 1
    tile2 = 2

    overlap_displacement(param_dir, tile1, tile2, plot=False)


if __name__ == "__main__":
    _main()
