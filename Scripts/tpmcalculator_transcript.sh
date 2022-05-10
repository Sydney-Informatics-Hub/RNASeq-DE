#!/bin/bash

sampleid=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
gtf=`echo $1 | cut -d ',' -f 3`
logfile=`echo $1 | cut -d ',' -f 4`
outdir=`echo $1 | cut -d ',' -f 5`

echo "$(date): Running TPMCalculator to obtain transript level TPM normalized counts. Sample ID:${sampleid}, BAM:${bam}, GTF: ${gtf}, Log file:${logfile}, Out:${outdir}" > ${logfile} 2>&1

cd ${outdir}

TPMCalculator -a -p -e \
	-q 255 \
	-g ${gtf} \
	-b ${bam} 2>>${logfile} 
