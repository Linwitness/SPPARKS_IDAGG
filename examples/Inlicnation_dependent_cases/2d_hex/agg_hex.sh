#!/bin/bash
fileBase="spparks"  #Sets base name for log, dump, cluster, and jpeg files
mpirun -np 1 ../../../src/spk_agg -var fileBase $fileBase < agg_hex.in
