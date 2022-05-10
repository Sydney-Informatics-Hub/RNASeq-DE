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

# Only ../${cohort}_final_bams/${sampleid}.final.bam are included
# using sampleids (column 2) from your cohort.config file

if [ -z "$1" ]
then
	echo "Please provide the path to your <cohort>.config file, e.g. sh read_distribution_make_input.sh cohort.config"
	exit
fi

config=$1
cohort=$(basename $config | cut -d'.' -f1)
logdir=./Logs/read_distribution
INPUTS=./Inputs
input_file=${INPUTS}/read_distribution.inputs
bamdir=../${cohort}_final_bams
outdir=../${cohort}_read_distribution

mkdir -p ${INPUTS} ${outdir} ${logdir}
rm -rf ${input_file}

tasks=()
while read -r fastq sampleid dataset reference seqcentre platform run_type library; do
	if [[ ! ${fastq} =~ ^#.*$ ]]; then
		bam=${bamdir}/${sampleid}.final.bam
		bed=$(ls ../Reference/${reference}/*${reference}*.bed)
		logfile=${logdir}/${sampleid}.log
		out=${outdir}/${sampleid}_read_distribution.txt
		tasks+=("${sampleid},${bam},${bed},${logfile},${out}")
	fi
done < "${config}"

unique=($(echo "${tasks[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
printf '%s\n' "${unique[@]}" > ${input_file}

num_tasks=`wc -l ${input_file}| cut -d' ' -f 1`
echo "Number of tasks: ${num_tasks}"
