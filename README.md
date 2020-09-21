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

Your directory structure should resemble the following: 

`.

├── Batch_1

│   ├── sample1_1.fastq.gz

│   └── sample1_2.fastq.gz

├── Batch_2

│   └── sample2.fastq.gz

├── README.md

├── References

│   ├── GRCh38

|   |   └── Homo_sapiens.GRCh38.dna.primary_assembly.fa

|   |   └── Homo_sapiens.GRCh38.94.gtf

│   └── GRCm38

|       └── Mus_musculus.GRCm38.dna.toplevel.fa

|       └── Mus_musculus.GRCm38.98.gtf

├── samples.config

└── Scripts

`

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

