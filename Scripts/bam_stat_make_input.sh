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

if [ -z "$1" ]
then
	echo "Please run this script with the path to your config file, e.g. sh bam_stat_make_input.sh cohort.config"
	exit
fi

config=$1
cohort=$(basename $config | cut -d'.' -f1)
logdir=./Logs/bam_stat
INPUTS=./Inputs
input_file=${INPUTS}/bam_stat.inputs
outdir=../${cohort}_final_bams

mkdir -p ${INPUTS} ${logdir}

rm -rf ${input_file}

tasks=()
while read -r fastq sampleid dataset reference seqcentre platform run_type library; do
	if [[ ! ${fastq} =~ ^#.*$ ]]; then
		bam=${outdir}/${sampleid}.final.bam
		logfile=${logdir}/${sampleid}.log
		out=${outdir}/${sampleid}_bam_stat.txt
		tasks+=("${sampleid},${bam},${logfile},${out}")
	fi
done < "${config}"

unique=($(echo "${tasks[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
printf '%s\n' "${unique[@]}" > ${input_file}

num_tasks=`wc -l ${input_file}| cut -d' ' -f 1`
echo "Number of tasks: ${num_tasks}"
