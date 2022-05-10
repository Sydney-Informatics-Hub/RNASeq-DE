#! /bin/bash

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

# All BAMs in ../${cohort}_final_bams/*bam are used
# This can include flowcell level BAMs

if [ -z "$1" ]
then
	echo "Please run this script <cohort>.config file, e.g. sh read_distribution_custom_make_input.sh <cohort>.config"
	exit
fi

config=$1
cohort=$(basename $config | cut -d'.' -f1)
bed=../Reference/GRCh38/Homo_sapiens.GRCh38.103.bed
logdir=./Logs/read_distribution
INPUTS=./Inputs
input_file=${INPUTS}/read_distribution.inputs
bamdir=../${cohort}_final_bams
outdir=../QC_reports/${cohort}_final_bams_read_distribution

mkdir -p ${INPUTS} ${logdir} ${outdir}

rm -rf ${input_file}

bams+=( $(ls $bamdir/*bam) )
for bam in "${bams[@]}"
do
	sampleid=$(basename $bam | cut -d'.' -f1)
	logfile=${logdir}/${sampleid}.log
	out=${outdir}/${sampleid}_read_distribution.txt
	echo "${sampleid},${bam},${bed},${logfile},${out}" >> ${input_file}
done

num_tasks=`wc -l ${input_file}| cut -d' ' -f 1`
echo "Number of tasks: ${num_tasks}"
