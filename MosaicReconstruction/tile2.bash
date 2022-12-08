echo "Block reconstruction for anatomix cerebellum tile #2"

echo -n "Started: "
date

time -v python block_reco_cerebellum_tile.py

echo -n "Finished: "
date

notify-send "Block reco tile #2 finished!"

mail -s "Block reco tile #2 finished!" mattia.humbel@unibas.ch << EOF
The reconstruction is done.
Cheers, Mattia
EOF
