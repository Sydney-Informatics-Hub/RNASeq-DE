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
	echo "Please run this script with the path to your config file, e.g. sh samtools_merge_index_make_input.sh ../cohort.config"
	exit
fi

config=$1
cohort=$(basename $config | cut -d'.' -f1)
logdir=./Logs/samtools_merge_index
outdir=../${cohort}_final_bams
INPUTS=./Inputs
input_file=${INPUTS}/samtools_merge_index.inputs

mkdir -p ${INPUTS} ${logdir} ${outdir}
rm -rf ${input_file}

# Ideally order tasks so multiplexed samples are merged first
# Then order by BAM size

samples=()
while read -r fastq sampleid dataset reference seqcentre platform run_type library; do
	if [[ ! ${fastq} =~ ^#.*$ ]]; then
		samples+=("${sampleid}")
	fi
done < "${config}"

unique=($(echo "${samples[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for sample in ${unique[@]}; do
	logfile=${logdir}/${sample}.log
	echo "${sample},${outdir},${logfile}" >> ${input_file}
done

num_tasks=`wc -l ${input_file}| cut -d' ' -f 1`
echo "Number of tasks: ${num_tasks}"
