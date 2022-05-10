#! /bin/bash

#########################################################
#
# Platform: NCI Gadi HPC
# Description: Run multiQC after running fastqc_run_parallel.pbs
# for <cohort>.config
# Usage: sh multiqc_all_datasets <cohort> 
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
#########################################################

# Summarizing fastQC report per unique dataset using multiQC
# Dataset name is obtained from column 3 of your cohort.config file

module load multiqc/1.9

if [ -z "$1" ]
then
	echo "Please run this script with the base name of your config file, e.g. sh multiqc_each_dataset.sh cohort.config"
	exit
fi

config=$1
cohort=$(basename "$config" | cut -d'.' -f 1)

# Get all unique datasets
while read -r fastq sampleid dataset reference seqcentre platform run_type library; do
	if [[ ! ${fastq} =~ ^#.*$ ]]; then
		fastqc=../${dataset}_fastQC
		datasets+=(${fastqc})	
	fi
done < "${config}"

# Unique datasets
uniq_datasets=($(echo "${datasets[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

echo "$(date): Performing multiqc for ${#uniq_datasets[@]} datasets: "${uniq_datasets[@]}""

for dataset in "${uniq_datasets[@]}"; do
	echo "$(date): MultiQC is running for ${dataset}..."
	multiqc --interactive ${dataset} -o ${dataset}
done
