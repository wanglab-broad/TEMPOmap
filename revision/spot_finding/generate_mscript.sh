#!/bin/bash
ppath="/stanley/WangLab/Connie/02.TEMPOmap/01.revision16Gene"
shpath="$ppath/code/spot_finding"
# listpath="$shpath/list"
outpath="$shpath/mscript"
main_mscript="$shpath/SF_2022_09_12_Rena_16Gene.m"

j=1
for i in {001..294}
# cat $listpath|while read file
do
	echo "tile='Position$i'" > $outpath/task_$j".m"
	cat $main_mscript >> $outpath/task_$j".m"
	j=$((j+1))
done
