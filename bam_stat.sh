#!/bin/bash

sampleid=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
logfile=`echo $1 | cut -d ',' -f 3`
out=`echo $1 | cut -d ',' -f 4`

module load python3/3.8.5
export PYTHONPATH=$HOME/.local/lib/python3.8/site-packages

echo "$(date): Running RSeQC's bam_stat.py to collect a summary of alignment metrics. Sample ID:${sampleid}, BAM:${bam}, Log file:${logfile}, Out:${out}" >> ${logfile} 2>&1

$HOME/.local/bin/bam_stat.py -i ${bam} > ${out} 2>${logfile} 
