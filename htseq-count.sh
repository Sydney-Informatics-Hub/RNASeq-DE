#!/bin/bash

sampleid=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
cohort=`echo $1 | cut -d ',' -f 3`
gtf=`echo $1 | cut -d ',' -f 4`
strand=`echo $1 | cut -d ',' -f 5`
logfile=`echo $1 | cut -d ',' -f 6`
NCPUS=`echo $1 | cut -d ',' -f 7`

module load python2

outdir=../${cohort}_htseq-count
out=${outdir}/${sampleid}.counts

mkdir -p ${outdir}
rm -rf ${logfile}

echo "$(date): Merging BAMs from multiplexed samples into a single BAM file. Sample ID:${sampleid}, Output:${outdir}, Log file: ${logfile}, NCPUS:${NCPUS}" >> ${logfile} 2>&1 

bam=${cohort}_final_bams/${sampleid}.final.bam

htseq-count -f bam -r pos --mode=union -s ${strand} ${bam} ${gtf} > ${out}
