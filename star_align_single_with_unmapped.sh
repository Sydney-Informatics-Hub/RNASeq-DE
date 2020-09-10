#!/bin/bash

# Align to the reference genome using STAR

filename=`echo $1 | cut -d ',' -f 1`
dataset=`echo $1 | cut -d ',' -f 2`
sampleid=`echo $1 | cut -d ',' -f 3`
fastq=`echo $1 | cut -d ',' -f 4`
ref=`echo $1 | cut -d ',' -f 5`
seqcentre=`echo $1 | cut -d ',' -f 6`
platform=`echo $1 | cut -d ',' -f 7`
library=`echo $1 | cut -d ',' -f 8`
lane=`echo $1 | cut -d ',' -f 9`
flowcell=`echo $1 | cut -d ',' -f 10`
outdir=`echo $1 | cut -d ',' -f 11`
logfile=`echo $1 | cut -d ',' -f 12`
NCPUS=`echo $1 | cut -d ',' -f 13`

echo `date` ": Mapping single reads with STAR 2.7.3a. Sample:$sampleid FASTQ:$fastq reference:$ref centre:$seqcentre platform:$platform lane:$lane flowcell:$flowcell NCPUS:$NCPUS" > ${logfile} 2>&1

# Mapping
STAR \
	--runThreadN ${NCPUS} \
	--genomeDir ${ref} \
	--quantMode GeneCounts \
	--readFilesCommand zcat \
	--readFilesIn ${fastq} \
	--outSAMattrRGline ID:${flowcell}:${lane} PU:${flowcell}.${lane}.${sampleid} SM:${sample} PL:${platform} CN:${seqcentre} LB:${library} \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx \
	--outFileNamePrefix ${outdir}/${sampleid}_${lane}_ >> ${logfile} 2>&1
