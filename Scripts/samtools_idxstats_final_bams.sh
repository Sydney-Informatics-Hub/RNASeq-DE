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

module load samtools/1.10
module load multiqc/1.9

if [ -z "$1" ]
then
        echo "Please provide the path to the directory containing BAMs, e.g. sh samtools_idxstats_final_bams.sh ../cohort_final_bams"
        exit
fi

bamdir=$(echo $1 | sed 's/\/$//')
outfileprefix=$(basename $bamdir)
echo $outfileprefix
outdir=../QC_reports/${outfileprefix}_samtools_idxstats
logdir=./Logs/samtools_idxstats

mkdir -p ${outdir} ${logdir}

# samtools_idxstats.sh in parallel on login node (48 parallel tasks)
find $bamdir -name "*bam" | xargs -i -n 1 -P 48 sh -c 'sample=$(basename {} | cut -d'.' -f1) && dir=$(basename $(dirname {} )) && samtools idxstats {} 1>../QC_reports/${dir}_samtools_idxstats/${sample}_idxstats.txt 2>./Logs/samtools_idxstats/$sample.log'

# Create multiqc report
multiqc --interactive -o ${outdir} ${outdir}
