#!/bin/bash

for subdir in */; do
    if [ -d "$subdir" ]; then
        echo "Processing $subdir"
        (cd "$subdir" && echo 1 13 | gmx hbond -s topol_prod.tpr -f traj_prod.xtc -n index.ndx -num hbond_counts.xvg -hbn hbonds.ndx -hbm hbonds.xpm -dt 1000)
    fi
done