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
readlen=`echo $1 | cut -d ',' -f 3`
logdir=`echo $1 | cut -d ',' -f 4`
adapters=`echo $1 | cut -d ',' -f 5`
NCPUS=`echo $1 | cut -d ',' -f 6`

basename=$(basename "$fastq" | cut -d. -f1)
logfile=${logdir}/${basename}_trimming.log

rm -rf ${logfile}

bbduk.sh -Xmx6g \
	threads=${NCPUS} \
	in=${fastq} \
	out=${out} \
	ref=${adapters} \
	ktrim=r \
	k=23 \
	mink=11 \
	hdist=1 \
	tpe \
	tbo \
	overwrite=true \
	trimpolya=${readlen} >> ${logfile} 2>&1
