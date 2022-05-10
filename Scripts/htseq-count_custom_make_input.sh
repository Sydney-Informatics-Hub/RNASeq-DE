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

# All bams here are included ../${cohort}_final_bams/*final.bam 

if [ -z "$1" ]
then
	echo "Please run this script your config file, e.g. sh htseq-count_custom_make_input.sh ../cohort.config"
	exit
fi

config=$1
cohort=$(basename $config | cut -d'.' -f1)
logdir=./Logs/htseq-count
strand=reverse
gtf=$(ls ../Reference/GRCh38/*gtf)
INPUTS=./Inputs
input_file=${INPUTS}/htseq-count.inputs

mkdir -p ${INPUTS} ${logdir}
rm -rf ${input_file}

tasks=()
bams+=( $(ls ../${cohort}_final_bams/*final.bam) )

for bam in "${bams[@]}"
do
	sampleid=$(basename $bam | cut -d'.' -f1)
	logfile=${logdir}/${sampleid}.log
	echo "${sampleid},${bam},${cohort},${gtf},${strand},${logfile}" >> $input_file
done

num_tasks=`wc -l ${input_file}| cut -d' ' -f 1`
echo "Number of tasks: ${num_tasks}"
