echo "Block reconstruction for anatomix cerebellum tile #2"

echo -n "Started: "
date

/usr/bin/time -v python block_reco_cerebellum_tile.py \
	/home/mattia/Documents/Cerebellum22/MosaicReconstruction/example/param_files/cerebellum_tile2.txt

echo -n "Finished: "
date

notify-send "Block reco tile #2 finished!"

mail -s "Block reco tile #2 finished!" mattia.humbel@unibas.ch << EOF
The reconstruction is done.
Cheers, Mattia
EOF
