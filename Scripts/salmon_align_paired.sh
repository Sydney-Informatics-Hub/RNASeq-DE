#!/bin/bash

#set -e

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

# Align paired reads to the reference transcriptome using salmon

filename=`echo $1 | cut -d ',' -f 1`
dataset=`echo $1 | cut -d ',' -f 2`
sampleid=`echo $1 | cut -d ',' -f 3`
fastq1=`echo $1 | cut -d ',' -f 4`
fastq2=`echo $1 | cut -d ',' -f 5`
ref=`echo $1 | cut -d ',' -f 6`
seqcentre=`echo $1 | cut -d ',' -f 7`
platform=`echo $1 | cut -d ',' -f 8`
library=`echo $1 | cut -d ',' -f 9`
lane=`echo $1 | cut -d ',' -f 10`
flowcell=`echo $1 | cut -d ',' -f 11`
outdir=`echo $1 | cut -d ',' -f 12`
logfile=`echo $1 | cut -d ',' -f 13`

#export PATH=$PATH:/scratch/er01/INFRA-121-CoreWorkflows-RNASeq/star-2.7.3a/bin/
#export PATH=$PATH:/g/data/er01/apps/star/2.7.3a/bin/

# Activate the conda environment
conda init bash
source ~/.bashrc
conda activate salmon


echo `date` ": Mapping with salmon. Sample:$sampleid R1:$fastq1 R2: $fastq2 centre:$seqcentre platform:$platform library:$library lane:$lane flowcell:$flowcell logfile:$logfile NCPUS:$NCPUS" > ${logfile} 2>&1

# Mapping
#STAR \
#        --runThreadN ${NCPUS} \
#        --outBAMsortingThreadN ${NCPUS} \
#        --genomeDir ${ref} \
#        --outBAMsortingBinsN 100 \
#        --quantMode GeneCounts \
#        --readFilesCommand zcat \
#        --readFilesIn ${fastq1} ${fastq2} \
#--outSAMattrRGline ID:${flowcell}:${lane} PU:${flowcell}.${lane}.${sampleid} SM:${sample} PL:${platform} CN:${seqcentre} LB:${library} \
#        --outSAMtype BAM SortedByCoordinate \
#        --outReadsUnmapped Fastx \
#        --outSAMunmapped Within KeepPairs \
#        --outFileNamePrefix ${outdir}/${sampleid}_${lane}_ >> ${logfile} 2>&1


#The <LIBTYPE> parameter should be replaced with the specific library type of the sequencing data. The library type describes the protocol used for generating the RNA-seq library and is important for accurate transcript quantification. Some common library types include:

# IU or ISR (unstranded or stranded): These library types do not indicate the strand of origin for the reads.
# ISF (stranded forward): The library is prepared in such a way that the sequencing reads originate from the forward (sense) strand of the transcripts.
# ISR (stranded reverse): The library is prepared in such a way that the sequencing reads originate from the reverse (antisense) strand of the transcripts.


# Nandan
	# need to re-place the salmon index in ../Referene/GRCh38 and remake the input file 
	# For now hard-coding it

# Quantifying in mapping-based mode
salmon quant -i /scratch/er01/INFRA-121-CoreWorkflows-RNASeq/salmon/salmon_index -l ISR -1 ${fastq1} -2 ${fastq2}  --validateMappings -o ${outdir}/${sampleid}_transcripts_quant -p 12 >> ${logfile} 2>&1



