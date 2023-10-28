#!/bin/bash
fileBase="spparks"  #Sets base name for log, dump, cluster, and jpeg files
ICfile="circleIC_000_000"
mpirun -np 1 ~/projects/SPPARKS-AGG/src/spk_agg -var fileBase $fileBase -var ICfile $ICfile < agg_embedded.in
