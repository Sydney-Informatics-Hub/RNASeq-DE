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

# Description: 
#	Creates input file for verify_md5_run_parallel.pbs
# 	for all samples in <cohort>.config. This is used to perform 
#  	checksums for files within or across datasets in parallel
#  	Expects checksum.md5 file to be in each dataset directory
#	1 task = 1 checksum for each line of checksum.md5
# Usage: 
# 	Create <cohort>.config file using template.
# 	sh verify_md5_make_input.sh <cohort>
# 	Then, use verify_md5_run_parallel.pbs
# 	Output: ../$dataset/verify_checksum_gadi.txt
# 	All files should have an OK, if not you need to investigate

if [ -z "$1" ]
then
	echo "Please provide the path to your config file, e.g. sh verify_md5_make_input.sh ../samples.config"
	exit
fi

config=$1
cohort=$(basename "$config" | cut -d'.' -f 1)
INPUTS=./Inputs
input_file=${INPUTS}/verify_md5.inputs

mkdir -p ${INPUTS}
rm -rf ${input_file}

datasets=($(awk '{print $3}' $config | uniq | xargs echo))
for dataset in "${datasets[@]}"; do
	#Ignore header
	if [[ ! $dataset =~ "DATASET" ]]; then
		checksum_file=../$dataset/checksums.md5
		# write to new output file
		rm -rf ../$dataset/verify_gadi_checksums.txt
		awk -v dataset=$dataset '{print dataset,$0}' $checksum_file >> ${input_file}
	fi
done

echo "Number of tasks: $(wc -l ${input_file} | awk '{print $1}')"
