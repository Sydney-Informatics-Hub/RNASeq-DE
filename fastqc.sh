#!/bin/bash

fastq=`echo $1 | cut -d ',' -f 1`
out=`echo $1 | cut -d ',' -f 2`
logfile=`echo $1 | cut -d ',' -f 3`
NCPUS=`echo $1 | cut -d ',' -f 4`

echo "$(date): Running fastQC. Fastq:${fastq}, Output:${out}, Log file: ${logfile}, NCPUS:${NCPUS}" >> ${logfile} 2>&1 

fastqc -t ${NCPUS} -o ${out} ${fastq} >> ${logfile} 2>&1 
