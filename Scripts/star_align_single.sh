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

# Align reads (not paired) to the reference genome using STAR
# `--outReadsUnmapped Fastx \` will retain unmapped reads

filename=`echo $1 | cut -d ',' -f 1`
dataset=`echo $1 | cut -d ',' -f 2`
sampleid=`echo $1 | cut -d ',' -f 3`
fastq=`echo $1 | cut -d ',' -f 4`
ref=`echo $1 | cut -d ',' -f 5`
seqcentre=`echo $1 | cut -d ',' -f 6`
platform=`echo $1 | cut -d ',' -f 7`
library=`echo $1 | cut -d ',' -f 8`
lane=`echo $1 | cut -d ',' -f 9`
flowcell=`echo $1 | cut -d ',' -f 10`
outdir=`echo $1 | cut -d ',' -f 11`
logfile=`echo $1 | cut -d ',' -f 12`

echo `date` ": Mapping single reads with STAR 2.7.3a. Sample:$sampleid FASTQ:$fastq reference:$ref centre:$seqcentre platform:$platform lane:$lane flowcell:$flowcell NCPUS:$NCPUS" > ${logfile} 2>&1

# Mapping
STAR \
	--runThreadN ${NCPUS} \
	--outBAMsortingThreadN ${NCPUS} \
	--genomeDir ${ref} \
	--quantMode GeneCounts \
        --outBAMsortingBinsN 100 \
	--readFilesCommand zcat \
	--readFilesIn ${fastq} \
	--outSAMattrRGline ID:${flowcell}:${lane} PU:${flowcell}.${lane}.${sampleid} SM:${sample} PL:${platform} CN:${seqcentre} LB:${library} \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx \
	--outFileNamePrefix ${outdir}/${sampleid}_${lane}_ >> ${logfile} 2>&1
