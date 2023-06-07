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

# Align paired reads to the reference genome using STAR
# Flags below will retain unmapped reads (as paired reads)
# 	--outReadsUnmapped Fastx \
#	--outSAMunmapped Within KeepPairs \

filename=`echo $1 | cut -d ',' -f 1`
dataset=`echo $1 | cut -d ',' -f 2`
sampleid=`echo $1 | cut -d ',' -f 3`
fastq1=`echo $1 | cut -d ',' -f 4`
fastq2=`echo $1 | cut -d ',' -f 5`
ref=`echo $1 | cut -d ',' -f 6`
seqcentre=`echo $1 | cut -d ',' -f 7`
platform=`echo $1 | cut -d ',' -f 8`
library=`echo $1 | cut -d ',' -f 9`
lane=`echo $1 | cut -d ',' -f 10`
flowcell=`echo $1 | cut -d ',' -f 11`
outdir=`echo $1 | cut -d ',' -f 12`
logfile=`echo $1 | cut -d ',' -f 13`

export PATH=$PATH:/scratch/er01/INFRA-121-CoreWorkflows-RNASeq/star-2.7.3a/bin/
#export PATH=$PATH:/g/data/er01/apps/star/2.7.3a/bin/


echo `date` ": Mapping with STAR. Sample:$sampleid R1:$fastq1 R2: $fastq2 centre:$seqcentre platform:$platform library:$library lane:$lane flowcell:$flowcell logfile:$logfile NCPUS:$NCPUS" > ${logfile} 2>&1

# Mapping
STAR \
        --runThreadN ${NCPUS} \
        --outBAMsortingThreadN ${NCPUS} \
        --genomeDir ${ref} \
        --outBAMsortingBinsN 100 \
        --quantMode GeneCounts \
        --readFilesCommand zcat \
        --readFilesIn ${fastq1} ${fastq2} \
        --outSAMattrRGline ID:${flowcell}:${lane} PU:${flowcell}.${lane}.${sampleid} SM:${sample} PL:${platform} CN:${seqcentre} LB:${library} \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outSAMunmapped Within KeepPairs \
        --outFileNamePrefix ${outdir}/${sampleid}_${lane}_ >> ${logfile} 2>&1
