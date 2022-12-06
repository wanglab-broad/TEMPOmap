#!/bin/bash -l
ppath="/stanley/WangLab/Connie/02.TEMPOmap/01.revision16Gene"
#$ -o /stanley/WangLab/Connie/02.TEMPOmap/01.revision16Gene/code/spot_finding/log/qsub_log_o.$JOB_ID.$TASK_ID
#$ -e /stanley/WangLab/Connie/02.TEMPOmap/01.revision16Gene/code/spot_finding/log/qsub_log_e.$JOB_ID.$TASK_ID

source "/broad/software/scripts/useuse"
reuse Matlab
now=$(date +"%T")
echo "Current time : $now"
command="run('$ppath/code/spot_finding/mscript/task_"$SGE_TASK_ID"');exit;"
matlab -nodisplay -nosplash -nodesktop -r $command

echo "Finished"
now=$(date +"%T")
echo "Current time : $now"
