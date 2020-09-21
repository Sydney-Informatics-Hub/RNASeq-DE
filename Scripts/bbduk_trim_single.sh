#!/bin/bash

fastq=`echo $1 | cut -d ',' -f 1`
out=`echo $1 | cut -d ',' -f 2`
readlen=`echo $1 | cut -d ',' -f 3`
logdir=`echo $1 | cut -d ',' -f 4`
adapters=`echo $1 | cut -d ',' -f 5`
NCPUS=`echo $1 | cut -d ',' -f 6`

basename=$(basename "$fastq" | cut -d. -f1)
logfile=${logdir}/${basename}_trimming.oe

rm -rf ${logfile}

echo "$(date): Trim 3' adapters and polyA tails from RNA seq - single read data" >> ${logfile} 2>&1
bbduk.sh -Xmx6g \
	threads=${NCPUS} \
	in=${fastq} \
	out=${out} \
	ref=${adapters} \
	ktrim=r \
	k=23 \
	mink=11 \
	hdist=1 \
	tpe \
	tbo \
	overwrite=true \
	trimpolya=${readlen} >> ${logfile} 2>&1

# These are the recommendations for adapter trimming by Brian Bushnell, one of the creators of BBDuk
# See http://seqanswers.com/forums/showthread.php?t=42776
# Added trimming of polya tail if they span the length of the entire read
# if your data is very low quality, you may wish to use more sensitive settings of hdist=2 k=21
# Can set memory by -Xmx3g, however when using the shell script, BBDuk can auto detect mem available
