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

fastq=`echo $1 | cut -d ',' -f 1`
out=`echo $1 | cut -d ',' -f 2`
logfile=`echo $1 | cut -d ',' -f 3`
NCPUS=`echo $1 | cut -d ',' -f 4`

echo "$(date): Running fastQC. Fastq:${fastq}, Output:${out}, Log file: ${logfile}, NCPUS:${NCPUS}" >> ${logfile} 2>&1 

fastqc -t ${NCPUS} -o ${out} ${fastq} >> ${logfile} 2>&1 
