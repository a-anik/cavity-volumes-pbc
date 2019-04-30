#!/bin/bash

rstart=0.1
dr=0.01
NSCAN=30
OUTFILE="Rscan_results.txt"

echo "# Rprobe Nvoids Vtotal Stotal" | tee "$OUTFILE"
for i in $(seq 0 $NSCAN);
do
    r=$(echo ${rstart}+${i}*${dr} | bc -l)
    cat Rscan_cavity_volumes_pbc.inp.tmpl | sed -e "s/@RPROBE@/$r/" > Rscan_cavity_volumes_pbc.inp
    ./cavity_volumes_pbc Rscan_cavity_volumes_pbc.inp
    N=`grep "Voids count" Rscan_cavout_info.dat  | cut -f2 -d":"`
    V=`grep "voids total volume" Rscan_cavout_info.dat  | cut -f2 -d":"`
    S=`grep "voids total surface" Rscan_cavout_info.dat  | cut -f2 -d":"`
    echo $r $N $V $S | tee -a "$OUTFILE"
done
