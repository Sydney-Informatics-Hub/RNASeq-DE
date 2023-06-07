#!/bin/bash

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

sampleid=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
cohort=`echo $1 | cut -d ',' -f 3`
gtf=`echo $1 | cut -d ',' -f 4`
strand=`echo $1 | cut -d ',' -f 5`
logfile=`echo $1 | cut -d ',' -f 6`
NCPUS=`echo $1 | cut -d ',' -f 7`

#module load python3/3.8.5
#module load python3/3.11.0


#export PYTHONPATH=$HOME/.local/lib/python3.8/site-packages
#export PYTHONPATH=/home/561/npd561/.local/lib/python3.8/site-packages
#export PYTHONPATH=/apps/python3/3.11.0/bin/python3.11/site-packages


#conda create -n htseq_env
conda init bash

# Activate the conda environment
source ~/.bashrc

conda activate htseq_env
#	conda init  bash
#	Restart the terminal
#htseq-count --version



echo $PYTHONPATH

outdir=../${cohort}_htseq-count
out=${outdir}/${sampleid}.counts

mkdir -p ${outdir}
rm -rf ${logfile}

echo "$(date): Running htseq-count to obtain raw counts. Sample ID:${sampleid}, BAM:${bam}, Cohort:${cohort}, Reference:${gtf}, Strand:${strand}, Output:${out}, Log file:${logfile}, NCPUS:${NCPUS}" >> ${logfile} 2>&1 

#$HOME/.local/bin/htseq-count -f bam -r pos --mode=union -s ${strand} ${bam} ${gtf} > ${out}

htseq-count -f bam -r pos --mode=union -s ${strand} ${bam} ${gtf} > ${out}


