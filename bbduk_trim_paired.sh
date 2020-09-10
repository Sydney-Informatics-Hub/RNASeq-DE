#!/bin/bash

module load bbtools/37.98

fastq1=`echo $1 | cut -d ',' -f 1`
fastq2=`echo $1 | cut -d ',' -f 2`
out1=`echo $1 | cut -d ',' -f 3`
out2=`echo $1 | cut -d ',' -f 4`
readlen=`echo $1 | cut -d ',' -f 5`
logdir=`echo $1 | cut -d ',' -f 6`
adapters=`echo $1 | cut -d ',' -f 7`
NCPUS=`echo $1 | cut -d ',' -f 8`

basename=$(basename "$fastq1" | cut -d. -f1)
uniq_basename="${basename::-1}"
logfile=${logdir}/${uniq_basename}trimming.oe

rm -rf ${logfile}

echo "$(date): Trim 3' adapters and polyA tails from RNA seq - paired read data" >> ${logfile} 2>&1

bbduk.sh -Xmx6g \
	threads=${NCPUS} \
	in=${fastq1} \
	in2=${fastq2} \
	out=${out1} \
	out2=${out2} \
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
