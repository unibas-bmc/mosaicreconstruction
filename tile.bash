#!/bin/bash

# Bash script to run the block reconstruction python program
# and send a notification, as well as an email, once it's done
# running.
#
# To be run, e.g. as:
# 	$ echo 'bash tile.bash' | batch
# or
#	$ nohup bash tile.bash &

tile_nr=4

echo "Block reconstruction for anatomix cerebellum tile #$tile_nr"

echo -n "Started: "
date

/usr/bin/time -v python block_reco_cerebellum_tile.py \
/home/mattia/Documents/Cerebellum22/MosaicReconstruction/\
example/param_files/cerebellum_tile$tile_nr.txt

echo -n "Finished: "
date

notify-send "Block reco tile #$tile_nr finished!"

mail -s "Block reco tile #$tile_nr finished!" mattia.humbel@unibas.ch \
<< EOF
The reconstruction is done.
Cheers, Mattia
EOF
