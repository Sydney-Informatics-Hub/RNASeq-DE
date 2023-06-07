#!/bin/bash

#set -e

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

#module load bbtools/37.98

fastq1=`echo $1 | cut -d ',' -f 1`
fastq2=`echo $1 | cut -d ',' -f 2`
out1=`echo $1 | cut -d ',' -f 3`
out2=`echo $1 | cut -d ',' -f 4`
readlen=`echo $1 | cut -d ',' -f 5`
logdir=`echo $1 | cut -d ',' -f 6`
adapters=`echo $1 | cut -d ',' -f 7`
NCPUS=`echo $1 | cut -d ',' -f 8`

basename=$(basename "$fastq1" | cut -d. -f1)
uniq_basename="${basename::-1}"
logfile=${logdir}/${uniq_basename}trimming.log


# Nandan
export PATH=$PATH:/scratch/er01/INFRA-121-CoreWorkflows-RNASeq/bbmap

rm -rf ${logfile}


#bbmap_PATH=/scratch/er01/INFRA-121-CoreWorkflows-RNASeq/bbmap

bbduk.sh -Xmx6g \
	threads=${NCPUS} \
	in=${fastq1} \
	in2=${fastq2} \
	out=${out1} \
	out2=${out2} \
	ref=${adapters} \
	ktrim=r \
	k=23 \
	mink=11 \
	hdist=1 \
	tpe \
	tbo \
	overwrite=true \
	trimpolya=${readlen} >> ${logfile} 2>&1

