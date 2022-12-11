#!/bin/bash
ppath="/stanley/WangLab/Connie/02.TEMPOmap/02.revisionCardiomyocyte64Gene"
shpath="$ppath/code/spot_finding"
# listpath="$shpath/list"
outpath="$shpath/mscript"
main_mscript="$shpath/spotFinding_2022_09_30_Rena_Cardiomyocyte64Gene.m"

j=1
for i in {001..528}
# cat $listpath|while read file
do
	echo "tile='Position$i'" > $outpath/task_$j".m"
	cat $main_mscript >> $outpath/task_$j".m"
	j=$((j+1))
done
