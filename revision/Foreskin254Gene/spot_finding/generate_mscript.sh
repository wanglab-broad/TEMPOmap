#!/bin/bash
ppath="/stanley/WangLab/Connie/02.TEMPOmap/04.revisionForeskin254Gene"
shpath="$ppath/code/spot_finding"
# listpath="$shpath/list"
outpath="$shpath/mscript"
main_mscript="$shpath/2022_10_15_Rena_SkinCulture254_gene.m"

j=1
for i in {001..432}
# cat $listpath|while read file
do
	echo "tile='Position$i'" > $outpath/task_$j".m"
	cat $main_mscript >> $outpath/task_$j".m"
	j=$((j+1))
done
