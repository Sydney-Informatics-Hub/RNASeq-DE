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

sampleid=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
gtf=`echo $1 | cut -d ',' -f 3`
logfile=`echo $1 | cut -d ',' -f 4`
outdir=`echo $1 | cut -d ',' -f 5`

echo "$(date): Running TPMCalculator to obtain TPM normalized counts. Sample ID:${sampleid}, BAM:${bam}, GTF: ${gtf}, Log file:${logfile}, Out:${outdir}" > ${logfile} 2>&1

cd ${outdir}

TPMCalculator -a -p -e \
	-q 255 \
	-g ${gtf} \
	-b ${bam} 2>>${logfile} 
