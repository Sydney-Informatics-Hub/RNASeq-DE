#! /bin/bash

# Create input file run gatk4 HaplotypeCaller in parallel
# Run before gatk4_hc_run_parallel.pbs

# Use bbtools provided list of adapters

module load bbtools/37.98

if [ -z "$1" ]
then
	echo "Please run this script with the base name of your config file, e.g. sh bbduk_trim_make_input.sh samples"
	exit
fi

cohort=$1
config=../$cohort.config
logs=./Logs/bbduk_trim
INPUTS=./Inputs
single=${INPUTS}/bbduk_trim_single.inputs
paired=${INPUTS}/bbduk_trim_paired.inputs
NCPUS=1

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
		

		echo "../${dataset}/$fastq,$out,$readlen,${logs},${adapters},${NCPUS}" >> ${single}
		
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
			echo "${fastq1},${fastq2},${out1},${out2},${readlen},${logs},${adapters},${NCPUS}" >> ${paired}
				
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



