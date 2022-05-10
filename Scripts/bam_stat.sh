#!/bin/bash

set -e

#########################################################
#
# Platform: NCI Gadi HPC
#
# Author: Tracy Chew
# tracy.chew@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the scientific and technical assistance
# <or e.g. bioinformatics assistance of <PERSON>> of Sydney Informatics
# Hub and resources and services from the National Computational
# Infrastructure (NCI), which is supported by the Australian Government
# with access facilitated by the University of Sydney.
#
#########################################################

module load python3/3.9.2
export PYTHONPATH=$HOME/.local/lib/python3.9/site-packages

sampleid=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
logfile=`echo $1 | cut -d ',' -f 3`
out=`echo $1 | cut -d ',' -f 4`

echo "$(date): Running RSeQC's bam_stat.py to collect a summary of alignment metrics. Sample ID:${sampleid}, BAM:${bam}, Log file:${logfile}, Out:${out}" >> ${logfile} 2>&1

$HOME/.local/bin/bam_stat.py -i ${bam} > ${out} 2>${logfile} 
