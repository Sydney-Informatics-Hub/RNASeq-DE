#!/bin/bash

sampleid=`echo $1 | cut -d ',' -f 1`
outdir=`echo $1 | cut -d ',' -f 2`
logfile=`echo $1 | cut -d ',' -f 3`
NCPUS=`echo $1 | cut -d ',' -f 4`

mkdir -p ${outdir}
rm -rf ${logfile}

echo "$(date): Merging BAMs from multiplexed samples into a single BAM file. Sample ID:${sampleid}, Output:${outdir}, Log file: ${logfile}, NCPUS:${NCPUS}" >> ${logfile} 2>&1 

sampbams=()
sampbams+=( $(ls ../${dataset}*_STAR/${sampleid}*.bam) )
final=${outdir}/${sampleid}.final.bam

if [ ${#sampbams[@]} -eq 0 ]; then
	echo "$(date): WARNING: No bams detected for ${sampleid}" >> ${logfile} 2>&1
elif [ ${#sampbams[@]} -eq 1 ];	then
	echo "$(date): ${sampleid} was not multiplexed. Renaming ${sampbams[0]} to ${final}" >> ${logfile} 2>&1
	cp ${sampbams[0]} ${final}
else
	echo "$(date): ${sampleid} was multiplexed. Merging ${sampbams[@]} bams into ${final}" >> ${logfile} 2>&1
	samtools merge -f -@ ${NCPUS} ${final} ${sampbams[@]} >> ${logfile} 2>&1
fi

echo "$(date): Indexing ${final}" >> ${logfile} 2>&1

samtools index -@ ${NCPUS} ${final} >> ${logfile} 2>&1
