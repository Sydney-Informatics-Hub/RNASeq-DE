#! /bin/bash

#########################################################
#
# Platform: NCI Gadi HPC
# Description: Creates input file for bam_stat_run_parallel.pbs 
# for all samples in <cohort>.config
# Usage: 
# sh bam_stat_make_inputs.sh <cohort>
# Author: Tracy Chew
# tracy.chew@sydney.edu.au
# Date last modified: 27/08/2020
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
########################################################

if [ -z "$1" ]
then
	echo "Please run this script with the base name of your ../<cohort>.config file, e.g. sh bam_stat_make_input.sh <cohort>"
	exit
fi

cohort=$1
config=../$cohort.config
logdir=./Logs/bam_stat
INPUTS=./Inputs
input_file=${INPUTS}/bam_stat.inputs
outdir=../${cohort}_final_bams

mkdir -p ${INPUTS}
mkdir -p ${logdir}

rm -rf ${input_file}

tasks=()
while read -r fastq sampleid dataset reference seqcentre platform run_type library; do
	if [[ ! ${fastq} =~ ^#.*$ ]]; then
		bam=${outdir}/${sampleid}.final.bam
		logfile=${logdir}/${sampleid}.oe
		out=${outdir}/${sampleid}_bam_stat.txt
		tasks+=("${sampleid},${bam},${logfile},${out}")
	fi
done < "${config}"

unique=($(echo "${tasks[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
printf '%s\n' "${unique[@]}" > ${input_file}

num_tasks=`wc -l ${input_file}| cut -d' ' -f 1`
echo "Number of samples to process for bam_stat_run_parallel.pbs: ${num_tasks}"
