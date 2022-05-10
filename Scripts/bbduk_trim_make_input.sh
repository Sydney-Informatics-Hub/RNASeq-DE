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

# Create inputs for bbduk_trim_run_parallel.pbs
# Will create separate input for single or paired reads (indicated in config file)
# Use bbtools provided list of adapters

if [ -z "$1" ]
then
        echo "Please provide the path to your config file, e.g. sh bbduk_trim_make_input.sh ../samples.config"
        exit
fi

config=$1
cohort=$(basename "$config" | cut -d'.' -f 1)
logs=./Logs/bbduk_trim
INPUTS=./Inputs
single=${INPUTS}/bbduk_trim_single.inputs
paired=${INPUTS}/bbduk_trim_paired.inputs

# set nocasematch option
shopt -s nocasematch

mkdir -p ${INPUTS}
mkdir -p ${logs}

rm -rf ${single}
rm -rf ${paired}

bbtools_path="$(module display bbtools | grep -m 1 PATH | awk '{print $3}')"
adapters=$bbtools_path/resources/adapters.fa

while read -r fastq sampleid dataset reference seqcentre platform run_type library; do
	# SINGLE end data
	if [[ ! ${fastq} =~ ^#.*$ && ${run_type} =~ single ]]; then
		outdir=../${dataset}_trimmed
		mkdir -p ${outdir}
		
		# Save trimmed fastq with same name as raw, except with _trimmed in prefix of file
		basename=$(basename "$fastq" | cut -d. -f1)
		single_extension="${fastq#*.}"
		out=${outdir}/${basename}_trimmed.${single_extension}
		
		# Get read lengths from 2nd line of each fastq file
		seq=`zcat "../${dataset}/${fastq}" | sed -n '2{p;q;}'`
		readlen=${#seq}
		
		singles+=(${fastq})
		

		echo "../${dataset}/$fastq,$out,$readlen,${logs},${adapters}" >> ${single}
		
	# PAIRED data
	elif [[ ! ${fastq} =~ ^#.*$ && ${run_type} =~ PAIRED ]]; then
		outdir=../${dataset}_trimmed
		mkdir -p ${outdir}
		
		# Only print to input file once per pair
		# Take basename of file, remove 1 or 2 from the end (typical naming convention of paired files)
		basename=$(basename "$fastq" | cut -d. -f1)
		paired_extension="${fastq#*.}"
		uniq_basename="${basename::-1}"
		which_pair="${basename: -1}"
		if [[ ${which_pair} =~ 1 ]]; then
			fastq1=../${dataset}/${uniq_basename}1.${paired_extension}
			fastq2=../${dataset}/${uniq_basename}2.${paired_extension}
			out1=${outdir}/${uniq_basename}1_trimmed.${paired_extension}
			out2=${outdir}/${uniq_basename}2_trimmed.${paired_extension}
			seq=`zcat "../${dataset}/${fastq1}" | sed -n '2{p;q;}'`
			readlen=${#seq}
			echo "${fastq1},${fastq2},${out1},${out2},${readlen},${logs},${adapters}" >> ${paired}
				
		fi
		pairs+=(${uniq_basename})
	fi
done < "${config}"

if [ ! -z "${singles}" ]
then
	echo "$(date): Found "${#singles[@]}" single read data in ${cohort}.config"
	echo "$(date): Saving single read input file for bbduk trimming to ${single}"
	
fi

if [ ! -z "${pairs}" ]
then
	uniq_pairs=($(printf "%s\n" "${pairs[@]}" | sort -u | tr '\n' ' '))
	echo "$(date): Found paired data. ASSUMES R1 and R2 IS DENOTED BY 1 AND 2 AT END OF FASTQ FILE PREFIX NAME"
	echo "$(date): Found "${#uniq_pairs[@]}" pairs in ${cohort}.config"
	echo "$(date): Saving paired data input file for bbduk trimming to ${paired}"
fi
