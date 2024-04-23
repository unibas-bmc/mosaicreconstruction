#!/bin/bash

basedir="/media/griffin/mouse4_perf_eth_db25/anatomix_mousebrain"
samplename="mouse4_perf_eth_db25"
testdir="$basedir/$samplename/stitchpos_tests"

workdir="/media/griffin/_home3/handover/MosaicReconstruction"
paramfile="$workdir/misc/elastix_rigid.txt"

fixed="$testdir/overlap1_reco.mhd"
moving="$testdir/overlap2_reco.mhd"
outdir="$testdir/overlap_elastix_out"
mask="$testdir/overlap_mask.mha"
mask_ratio="0.6"

./createCircMask.py $fixed $mask $mask_ratio

if [[ ! -d $outdir ]]
then
	mkdir $outdir
fi

elastix -f $fixed -m $moving -fMask $mask -mMask $mask -out $outdir \
       	-p $paramfile

echo "(TransformParameters theta_x, theta_y, theta_z, t_x, t_y, tz)"
echo "angles in radian, translation in pixels"
grep "(TransformParameters " $outdir/TransformParameters.0.txt
