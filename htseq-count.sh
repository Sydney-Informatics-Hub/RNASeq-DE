#!/bin/bash

sampleid=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
cohort=`echo $1 | cut -d ',' -f 3`
gtf=`echo $1 | cut -d ',' -f 4`
strand=`echo $1 | cut -d ',' -f 5`
logfile=`echo $1 | cut -d ',' -f 6`
NCPUS=`echo $1 | cut -d ',' -f 7`

module load python3/3.8.5
export PYTHONPATH=$HOME/.local/lib/python3.8/site-packages

echo $PYTHONPATH

outdir=../${cohort}_htseq-count
out=${outdir}/${sampleid}.counts

mkdir -p ${outdir}
rm -rf ${logfile}

echo "$(date): Running htseq-count to obtain raw counts. Sample ID:${sampleid}, BAM:${bam}, Cohort:${cohort}, Reference:${gtf}, Strand:${strand}, Output:${out}, Log file:${logfile}, NCPUS:${NCPUS}" >> ${logfile} 2>&1 

$HOME/.local/bin/htseq-count -f bam -r pos --mode=union -s ${strand} ${bam} ${gtf} > ${out}
