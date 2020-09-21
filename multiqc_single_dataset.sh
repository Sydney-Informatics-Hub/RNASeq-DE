#! /bin/bash

module load python2

if [ -z "$1" ]
then
	echo "Please run this script with the base name of your ../<cohort>.config file, e.g. sh multiqc.sh <config>"
	exit
fi

cohort=$1
fastqc=../${cohort}\_fastQC

echo "Running multiQC for fastQC files in ${fastqc}..."
multiqc --interactive ${fastqc} -o ${fastqc}
