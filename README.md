# RNASeq-DE

The scripts in this repository process RNA sequencing data (single, paired and/or multiplexed) for differential expression (raw FASTQ to counts) on the National Compute Infrastructure, Gadi. The scripts include:

1. FastQC to obtain quality reports on fastq files
2. MultiQC to summaries quality reports on fastq files
3. BBduk trim to trim 3' adapters and poly A tails from raw FASTQ files
4. STAR for spliced aware alignment of RNA sequencing reads to a reference genome
6. SAMtools to merge and index sample BAMs
5. RSeQC for obtaining a summary of alignment metrics from BAM files
6. HTSeq for obtaining raw counts 

# Set up

The scripts in this repository use relative paths and require careful setup to ensure that scripts can locate input files seamlessly. To start, go to your working directory (e.g. /scratch/ab1 for your ab1 project on NCI Gadi) and:

`git clone https://github.com/Sydney-Informatics-Hub/RNASeq-DE`

`cd RNASeq-DE`

You will then need to include:

* Your __Raw FASTQs__ into this directory, keeping FASTQ files from 1 sequencing batch in 1 unique "batch" directory
* Your __reference files__. Reference genome primary assembly (.fasta) and corresponding annotation (.gtf) file needs to be in the relevant `Reference` sub-directory. References can be obtained:
    * following recommendations in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). 
    * from [Ensembl](https://asia.ensembl.org/info/data/ftp/index.html) 
    * SIH's CloudStor, which contains already genomes already indexed for STAR v2.7.3a (coming soon)!

Your __RNASeq-DE__ directory structure should resemble the following: 

```bash
├── Batch_1
│   ├── sample1_1.fastq.gz
│   └── sample1_2.fastq.gz
├── Batch_2
│   └── sample2.fastq.gz
├── README.md
├── References
│   ├── GRCh38
│   │   └── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   │   └── Homo_sapiens.GRCh38.94.gtf
│   └── GRCm38
│       └── Mus_musculus.GRCm38.dna.toplevel.fa
│       └── Mus_musculus.GRCm38.98.gtf
├── samples.config
└── Scripts
```

Edit the __samples.config__ file. An example is provided below:

|#FASTQ|	SAMPLEID|	DATASET|	REFERENCE_GRCh38_GRCm38|	SEQUENCING_CENTRE|	PLATFORM|	RUN_TYPE_SINGLE_PAIRED|	LIBRARY|
|------|---------|----------|------------------------|--------------------|-----------|-------------------------|--------|
|sample1_1.fastq.gz|	SAMPLEID1|	Batch_1|	GRCh38|	KCCG|	ILLUMINA|	PAIRED|	1|
|sample1_2.fastq.gz|     SAMPLEID1|       Batch_1| GRCh38|  KCCG|    ILLUMINA|        PAIRED  1|
|sample2.fastq.gz|	SAMPLEID2|	Batch_2|	GRCm38|	KCCG|	ILLUMINA|	SINGLE|	1|

Column descriptions for __samples.config__:

|Column name| Description|
|----|--------|
|FASTQ| FASTQ file name. This column can be populated with `ls -1` in your sequencing batch directory|
|SAMPLEID| The sample identifier used in your laboratory. This will be used in naming output files. No whitespace please.|
|DATASET| The sequencing batch that the FASTQ file was generated, and the directory name where the FASTQ file is located. No whitespace please. |
|REFERENCE_GRCh38_GRCm38| Reference subdirectory name, e.g. GRCh38 or GRCm38 in the above example. Scripts will use reference files (.fasta and .gtf) and STAR index files for the FASTQ file/sample for alignment and counting.  |
|SEQUENCING_CENTRE| e.g. KCCG. This is used in the read group header in the output BAM file for the aligned FASTQ. No whitespace please.|
|PLATFORM| e.g. ILLUMINA. This is used in the read group header in the output BAM file for the aligned FASTQ.|
|RUN_TYPE_SINGLE_PAIRED| Whether you want to trim, STAR align the FASTQ as single read data or paired end data.| 
|LIBRARY| The sequencing library of the FASTQ file. This is used in the read group header in the output BAM file for the aligned FASTQ. No whitespace please.|

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

