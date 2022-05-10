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

module load multiqc/1.9

if [ -z "$1" ]
then
	echo "Please run this script with the path to input directory e.g. sh multiqc.sh ../dataset_fastQC"
	exit
fi

dir=$1

echo "Running multiQC for files in ${dir}..."
multiqc --interactive ${dir} -o ${dir}
