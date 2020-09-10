#!/bin/bash

# Inputs required by star single: reference, sampleid, fastq, out
# Inputs requred by star paired: reference, sampleid, fastq1, fastq2, out
# Assumes trimming was performed with bbduk_trim_run_parallel.pbs
# Need to sort to align fastq with most reads - least reads for CPU efficiency

if [ -z "$1" ]
then
	echo "Please run this script with the base name of your config file, e.g. sh star_align_trimmed_make_input.sh samples"
	exit
fi

cohort=$1
config=../$cohort.config
INPUTS=./Inputs
single=${INPUTS}/star_align_single.inputs
paired=${INPUTS}/star_align_paired.inputs
unsort_single=${INPUTS}/star_single_unsort.inputs
unsort_paired=${INPUTS}/star_paired_unsort.inputs
singletmp=${INPUTS}/star_single.tmp
pairedtmp=${INPUTS}/star_paired.tmp
NCPUS=12	
# 24 CPUs allows 96 Gb mem on Gadi normal nodes
# 4 CPUs allows 125Gb mem on Gadi hugemem nodes

mkdir -p ${INPUTS}

rm -rf ${single}
rm -rf ${paired}
rm -rf ${unsort_single}
rm -rf ${unsort_paired}
rm -rf ${singletmp}
rm -rf ${pairedtmp}

# set nocasematch option
shopt -s nocasematch

while read -r fastq sampleid dataset reference seqcentre pl run_type library; do
	ref=../Reference/${reference}
	# Convert platform name to all uppercase
	platform=${pl^^}
	
	if [[ ! ${fastq} =~ ^#.*$ && ${run_type} =~ single ]]; then
		indir=../${dataset}_trimmed
		outdir=../${dataset}_STAR
		logdir=./Logs/star_align_trimmed/${dataset}
			
		mkdir -p ${outdir}
		mkdir -p ${logdir}

		# Save trimmed fastq with same name as raw, except with _trimmed in prefix of file
		basename=$(basename "$fastq" | cut -d. -f1)	# use basename to sort files by read number
		single_extension="${fastq#*.}"
		fastq1=${indir}/${basename}_trimmed.${single_extension}
		
		# Get flowcell and lane from fastq, assuming typical Illumina naming convention (Illumina pipelines 1.4 onwards)
		#flowcell=$(zcat ${fastq1} | head -1 | cut -d ':' -f 3)
		#lane=$(zcat ${fastq1} | head -1 | cut -d ':' -f 4)

		# UBC sequenced samples have weird convention...there is no flowcell ID
		flowcell=NA
                lane=$(zcat ${fastq1} | head -1 | cut -d ':' -f 2)

		# Log file to save to
		logfile=${logdir}/${sampleid}_${lane}.oe
				
		echo "${basename},${dataset},${sampleid},${fastq1},${ref},${seqcentre},${platform},${library},${lane},${flowcell},${outdir},${logfile},${NCPUS}" >> ${unsort_single}
		
		num_single+=("${sampleid}")
		num_samples+=("${sampleid}")
		single_datasets+=("${dataset}")
		
		
	elif [[ ! ${fastq} =~ ^#.*$ && ${run_type} =~ PAIRED ]]; then
		indir=../${dataset}_trimmed
		outdir=../${dataset}_STAR
		logdir=./Logs/star_align_trimmed/${dataset}

		mkdir -p ${outdir}
		mkdir -p ${logdir}
		
		basename=$(basename "$fastq" | cut -d. -f1)
		paired_extension="${fastq#*.}"
		uniq_basename="${basename::-1}"
		which_pair="${basename: -1}"
		
		if [[ ${which_pair} =~ 1 ]]; then
			fastq1=${indir}/${uniq_basename}1_trimmed.${paired_extension}
			fastq2=${indir}/${uniq_basename}2_trimmed.${paired_extension}
			#seq=`zcat "../${dataset}/${fastq1}" | sed -n '2{p;q;}'`
			
			# Get flowcell and lane from fastq, assuming typical illumina naming convention
			flowcell=$(zcat ${fastq1} | head -1 | cut -d ':' -f 3)
			lane=$(zcat ${fastq1} | head -1 | cut -d ':' -f 4)
			
			# Log file to save to
			logfile=${logdir}/${sampleid}_${lane}.oe
			
			echo "${basename},${dataset},${sampleid},${fastq1},${fastq2},${ref},${seqcentre},${platform},${library},${lane},${flowcell},${outdir},${logfile},${NCPUS}" >> ${unsort_paired}
			num_paired+=("${sampleid}")
		fi
		num_samples+=("${sampleid}")
		paired_datasets+=("${dataset}")		
	fi
done < "${config}"


echo "$(date): $config has ${#num_samples[@]} samples. Number of samples with single read data: ${#num_single[@]}. Number of samples with paired read data: ${#num_paired[@]}."

# Get number of reads per single read fastq file and total reads per fastq pair 
# Order from highest to lowest
# Use multiqc_data/multiqc_general_stats.txt
uniq_single_datasets=($(printf "%s\n" "${single_datasets[@]}" | sort -u | tr '\n' ' '))
uniq_paired_datasets=($(printf "%s\n" "${paired_datasets[@]}" | sort -u | tr '\n' ' '))

# Sort single input by number of reads 
if (( ${#uniq_single_datasets[@]} )); then
	for dataset in "${uniq_single_datasets[@]}"
	do
	        echo "${dataset} has single read data"
        	multiqc=../${dataset}_fastQC/multiqc_data/multiqc_general_stats.txt
	        while read -r sampleid dups readlen total fail gc; do
        	if [[ ! $sampleid =~ Sample ]]; then
	                # Warning - sampleid wasn't unique in this dataset -_-
                	# Match by dataset or sequencing batch too
        	        echo "$sampleid,$dataset,$total" >> ${singletmp}
	        fi
        	done < "${multiqc}"
	done
	sort -k2 -t, -rg ${singletmp} -o ${singletmp}

	# use $sample and $dataset in $singletmp to sort $single
	# need both because $sample = filename used by seq company which is not always unique
	for S in $(cat ${singletmp} | awk -F, '{print $1","$2}'); do grep -P ^$S, ${unsort_single} >> ${single}; done
fi 

if (( ${#uniq_paired_datasets[@]} )); then
	# Sort paired input by number of reads 
	for dataset in "${uniq_paired_datasets[@]}"
	do
		echo "${dataset} has paired data"
		multiqc=../${dataset}_fastQC/multiqc_data/multiqc_general_stats.txt
		while read -r sampleid dups readlen total fail gc; do
		if [[ $sampleid =~ 1$ && ! $sampleid =~ Sample ]]; then
			total_per_pair=$(expr $total*2 | bc)
			echo "$sampleid,$dataset,$total_per_pair" >> ${pairedtmp}
		fi
		done < "${multiqc}"
	done
	sort -k2 -t, -rg ${pairedtmp} -o ${pairedtmp}
	for S in $(cat ${pairedtmp} | awk -F, '{print $1","$2}'); do grep -P ^$S, ${unsort_paired} >> ${paired}; done
fi

rm -rf ${unsort_single}
rm -rf ${unsort_paired}
rm -rf ${singletmp}
rm -rf ${pairedtmp}
