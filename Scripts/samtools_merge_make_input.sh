#! /bin/bash

#########################################################
#
# Platform: NCI Gadi HPC
# Description: Creates input file for samtools_merge_run_parallel.pbs 
# for all samples in <cohort>.config
# Usage: 
# sh samtools_merge_make_input.sh <cohort>
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
	echo "Please run this script with the base name of your config file, e.g. sh samtools_merge_make_input.sh samples"
	exit
fi

cohort=$1
config=../$cohort.config
logdir=./Logs/samtools_merge
outdir=../${cohort}_final_bams
INPUTS=./Inputs
input_file=${INPUTS}/samtools_merge.inputs
NCPUS=1

mkdir -p ${INPUTS}
mkdir -p ${logdir}
mkdir -p ${outdir}

rm -rf ${input_file}

samples=()
while read -r fastq sampleid dataset reference seqcentre platform run_type library; do
	if [[ ! ${fastq} =~ ^#.*$ ]]; then
		samples+=("${sampleid}")
	fi
done < "${config}"

unique=($(echo "${samples[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for sample in ${unique[@]}; do
	logfile=${logdir}/${sample}.oe
	echo "${sample},${outdir},${logfile},${NCPUS}" >> ${input_file}
done

num_tasks=`wc -l ${input_file}| cut -d' ' -f 1`
echo "Number of samples to process for samtools_merge_run_parallel.pbs: ${num_tasks}"
