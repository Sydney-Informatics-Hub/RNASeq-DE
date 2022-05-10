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
outdir=`echo $1 | cut -d ',' -f 2`
logfile=`echo $1 | cut -d ',' -f 3`

mkdir -p ${outdir}
rm -rf ${logfile}

sampbams=()
sampbams+=( $(ls ../*STAR/${sampleid}_*.bam) )
final=${outdir}/${sampleid}.final.bam

if [ ${#sampbams[@]} -eq 0 ]; then
	echo "$(date): WARNING: No bams detected for ${sampleid}" >> ${logfile} 2>&1
elif [ ${#sampbams[@]} -eq 1 ];	then
	echo "$(date): ${sampleid} was not multiplexed. Renaming ${sampbams[0]} to ${final}" >> ${logfile} 2>&1
	ln -sf ${sampbams[0]} ${final}
	echo "$(date): Indexing ${final}" >> ${logfile} 2>&1
	samtools index -@ ${NCPUS} ${final} >> ${logfile} 2>&1
else
	# More than 2 bams for sample found. Merge first and index merged BAM first (takes longer)
        echo "$(date): ${sampleid} was multiplexed. Merging ${sampbams[@]} bams into ${final}" >> ${logfile} 2>&1
        samtools merge -f -@ ${NCPUS} ${final} ${sampbams[@]} >> ${logfile} 2>&1
        echo "$(date): Indexing ${final}" >> ${logfile} 2>&1
        samtools index -@ ${NCPUS} ${final} >> ${logfile} 2>&1

        echo "$(date): Retaining per batch BAMs for ${sampleid}" >> ${logfile} 2>&1
        for batchbam in "${sampbams[@]}"
        do
                echo "$(date): Renaming ${batchbam} to ${outdir}/${sampleid}_${dataset}.final.bam" >> ${logfile} 2>&1
                dataset=`echo $batchbam | cut -d'/' -f2 | cut -d'_' -f 3` ## CUSTOM FOR PREDICT-19
                # Deleting symlink created with ln also deletes original. Much safer to use cp -s method which keeps original
		ln -sf ${batchbam} ${outdir}/${sampleid}_${dataset}.final.bam
                # Can't use relative paths...
		#cp -s ${batchbam} ${outdir}/${sampleid}_${dataset}.final.bam
		echo "$(date): Indexing ${outdir}/${sampleid}_${dataset}.final.bam" >> ${logfile} 2>&1
                samtools index -@ ${NCPUS} ${outdir}/${sampleid}_${dataset}.final.bam >> ${logfile} 2>&1
        done
fi
