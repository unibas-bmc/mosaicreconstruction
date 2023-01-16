#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import norm

class RecoParam:
	def __init__(paramfile):
        # read file, make structure array
        this.infoStruct = {}
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
        this.samplename = infoStruct['samplename']
        this.infofile = infoStruct['infofile']
        this.projbasedir = infoStruct['projpath2']
        this.recobasedir = infoStruct['recopath']
        this.stripheight = float(infoStruct['stripheight'])
        this.outputcrop = infoStruct['outputcrop'].split(",")
        for i in range(len(this.outputcrop)):
            this.outputcrop[i] = int(this.outputcrop[i])

        this.outputgrayscale = infoStruct['outputgrayscale'].split(",")
        for i in range(len(this.outputgrayscale)):
            this.outputgrayscale[i] = float(this.outputgrayscale[i])

        this.blocksize = int(infoStruct['stripheight'])
        this.ncore = int(infoStruct['ncore'])
        this.pixsize_um = float(infoStruct['pixsize_um'])


param1 = RecoParam("/home/mattia/Documents/Cerebellum22/MosaicReconstruction/example/param_files/cerebellum_tile2.txt")
param2 = RecoParam("/home/mattia/Documents/Cerebellum22/MosaicReconstruction/example/param_files/cerebellum_tile3.txt")


width = 14000. # pixels
radius = 7500. # pixels

pixsize = 6.5e-4 # mm

# u & v motor positions
x1 = np.array([0.03000608, 8.26400422]) # mm
x2 = np.array([7.22999703, 4.107003875]) # mm

x1 /= pixsize
x2 /= pixsize

x2 -= x1
x1 = np.array([np.ceil(width / 2), np.ceil(width/2)])
x2 += x1

center = np.round((x1 + x2) / 2.)
rect = np.broadcast_to(center, (4, 2))

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

while rect_is_contained(tmp := rect_grow(rect), x1, x2, radius):
	rect = tmp

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

rect2 = np.round(rect - (x2 - x1).reshape(1,2))
# assume order z, y, x
s1 = np.s_[rect[0,1]:rect[3,1], rect[0,0]:rect[1,0]]
s2 = np.s_[rect2[0,1]:rect2[3,1], rect2[0,0]:rect2[1,0]]

print(s1)
print(s2)

# plt.show()
