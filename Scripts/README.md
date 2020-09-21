# RNASeq-DE

The scripts in this repository contain scripts that process RNA sequencing data (single, paired and/or multiplexed) for differential expression (raw FASTQ to counts) on the National Compute Infrastructure, Gadi. The scripts include:

1. FastQC to obtain quality reports on fastq files
2. MultiQC to summaries quality reports on fastq files
3. BBduk trim to trim 3' adapters and poly A tails from raw FASTQ files
4. STAR for spliced aware alignment of RNA sequencing reads to a reference genome
6. SAMtools to merge and index sample BAMs
5. RSeQC for obtaining a summary of alignment metrics from BAM files
6. HTSeq for obtaining raw counts 

# Set up

## References

Download your reference genome primary assembly and corresponding annotation file (GTF format) following recommendations in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). You can download the FASTA and GTF file from from an online repository such as [Ensembl](https://asia.ensembl.org/info/data/ftp/index.html) for your species. 

## Software

openmpi/4.0.2
nci-parallel/1.0.0
fastQC/0.11.7
python3/3.8.5
`pip3 install multiqc`
BBDuk/37.98
STAR/2.7.3a
SAMtools/1.10
`pip3 install RSeQC`
`pip3 install HTSeq`
