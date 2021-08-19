# RNASeq-DE

The scripts in this repository process RNA sequencing data (single, paired and/or multiplexed) for differential expression (raw FASTQ to counts) on the __National Compute Infrastructure, Gadi__. The scripts include:

1. FastQC to obtain quality reports on raw fastq files
2. MultiQC to summarize fastQC quality reports on raw fastq files
3. BBduk trim to trim 3' adapters and poly A tails from raw FASTQ files
4. FastQC to obtain quality reports on trimmed fastq files
5. MultiQC to summarize fastQC quality reports on trimmed fastq files
6. STAR indexing for a given reference genome and corresponding annotation file
7. STAR for spliced aware alignment of trimmed RNA sequencing reads to a reference genome
8. SAMtools to merge and index sample BAMs
9. RSeQC for obtaining a summary of alignment metrics from BAM files
10. HTSeq for obtaining raw counts 

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

Edit the __cohort.config__ file. This config file is used to tell the scripts which samples to process, how to process them, and where it can locate relevant input files. An example is provided below:

|#FASTQ|	SAMPLEID|	DATASET|	REFERENCE_GRCh38_GRCm38|	SEQUENCING_CENTRE|	PLATFORM|	RUN_TYPE_SINGLE_PAIRED|	LIBRARY|
|------|---------|----------|------------------------|--------------------|-----------|-------------------------|--------|
|sample1_1.fastq.gz|	SAMPLEID1|	Batch_1|	GRCh38|	KCCG|	ILLUMINA|	PAIRED|	1|
|sample1_2.fastq.gz|     SAMPLEID1|       Batch_1| GRCh38|  KCCG|    ILLUMINA|        PAIRED|  1|
|sample2.fastq.gz|	SAMPLEID2|	Batch_2|	GRCm38|	KCCG|	ILLUMINA|	SINGLE|	1|

Column descriptions for __cohort.config__:

|Column name| Description|
|----|--------|
|FASTQ| FASTQ file name. This column can be populated with `ls -1` in your sequencing batch directory|
|SAMPLEID| The sample identifier used in your laboratory. This will be used in naming output files. Avoid whitespace.|
|DATASET| The sequencing batch that the FASTQ file was generated, and the directory name (must case match) where the FASTQ file is located. |
|REFERENCE_GRCh38_GRCm38| Reference subdirectory name, e.g. GRCh38 or GRCm38 in the above example (must case match). Scripts will use reference files (.fasta and .gtf) and STAR index files for the FASTQ file/sample for alignment and counting.  |
|SEQUENCING_CENTRE| e.g. KCCG. This is used in the read group header in the output BAM file for the aligned FASTQ. No whitespace please.|
|PLATFORM| e.g. ILLUMINA. This is used in the read group header in the output BAM file for the aligned FASTQ.|
|RUN_TYPE_SINGLE_PAIRED| Input SINGLE or PAIRED. This is used to indicate whether you want to trim, STAR align the FASTQ as single read data or paired end data.| 
|LIBRARY| The sequencing library of the FASTQ file. This is used in the read group header in the output BAM file for the aligned FASTQ. No whitespace please.|

## Software

The software listed below are used in the RNASeq-DE pipeline. Some of these are installed globally on NCI Gadi (check with `module avail` for the current software). Install python3 packages by `module load python3/3.8.5`, and then using the `pip3 install` commands as listed below. These will be installed in `$HOME`. All other software need to be installed in your project's `/scratch` directory and module loadable. 

openmpi/4.0.2 (installed globally)

nci-parallel/1.0.0 (installed globally)

SAMtools/1.10 (installed globally)

python3/3.8.5 (installed globally)

`pip3 install multiqc`

`pip3 install RSeQC`

`pip3 install HTSeq`

fastQC/0.11.7

BBDuk/37.98

STAR/2.7.3a

# Running the pipeline

Change into the Scripts directory using `cd Scripts` and run all steps below in sequence to process data for samples belonging to `dataset` in `..cohort.config`. 

Generally, each step will involve running a make_input script first. This makes inputs for the run_parallel.pbs script which is run next, once you have modified compute resource requests according to your data requirements (a guide is provided in the script and down below). The run_parallel.pbs script runs a task.sh script in parallel with the requested compute resources per job and per task, where 1 row in the input file = 1 task. 

1. Obtain quality reports with FastQC for your raw fastq files by:
      * `sh fastqc_make_input.sh cohort`. 
      * Editing `fastqc_run_parallel.pbs` project id and compute requirements and submitting the job by `qsub fastqc_run_parallel.pbs`
      * Description: `fastqc_make_input.sh` creates inputs for each fastq file. FastQC reports are written to `../dataset_fastQC` (1 fastq file = 1 `fastqc.sh` task). Inputs are ordered by filesize.
2. Summarize fastQC quality reports for raw fastq files using MultiQC by:
      * `sh multiqc.sh cohort`
      * Description: This will run multiqc for each raw fastq for dataset in `..cohort.config`. This job only takes a couple of seconds (for cohorts with <100 samples, ~30M paired end reads) and can be run on the command line. MultiQC reports are written to `../dataset_fastQC`. Before moving on, you should assess your multiqc output by looking at `multiqc_report.html` on a browser or checking `multiqc_data/multiqc_fastqc.txt and multiqc_data/multiqc_general_stats.txt`
3. Trim 3' adapters and poly A tails from raw FASTQ files by:
      * `sh bbduk_trim_make_input.sh cohort`
      * Editing `bbduk_trim_run_parallel.pbs` project id and compute requirements and submitting the job by `qsub bbduk_trim_run_parallel.pbs`
      * Description: This script performs trimming of raw FASTQ files for single read and/or paired read data as specified by the RUN_TYPE_SINGLE_PAIRED column in your config file. Trimmed FASTQs are written to `../dataset_trimmed`. By default, scripts perform adapter trimming following recommendations under "Adapter Trimming" on the [BBDuk Guide](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/), and in addition, will perform poly A tail trimming for the entire length of the input reads.
4. Obtain quality reports with FastQC for trimmed fastq files by:
      * `sh fastqc_trimmed_make_input.sh cohort`
      * Editing `fastqc_trimmed_run_parallel.pbs` project id and compute requirements and submitting the job by `qsub fastqc_trimmed_run_parallel.pbs`. 
      * Description: `fastqc_trimmed_make_input.sh` creates inputs for each trimmed fastq file in `../dataset_trimmed`. FastQC reports are written to `../dataset_fastQC_trimmed` (1 fastq file = 1 `fastqc.sh` task). Inputs are ordered by filesize.
 5. Summarize fastQC quality reports for trimmed fastq files using MultiQC by:
      * `sh multiqc_trimmed.sh cohort`
      * This will run multiqc for each trimmed fastq file for dataset in `..cohort.config`. This job only takes a couple of seconds (for cohorts with <100 samples, ~30M paired end reads) and can be run on the command line. MultiQC reports are written to `../dataset_fastQC_trimmed`. Before moving on, you should assess your multiqc output by looking at `multiqc_report.html` on a browser or checking `multiqc_data/multiqc_fastqc.txt and multiqc_data/multiqc_general_stats.txt` and compare it to quality reports for the raw FASTQ files.
6. Index reference genome and annotation file using STAR. 
      * Description: Indexing must be performed using the same version of STAR that will be used for the STAR mapping step
      
# Benchmarking metrics

The below metrics were obtained for a mouse dataset with 10 samples, ~33 M reads, 150 base pair, paired end reads. 

|#JobName|CPUs_requested|CPUs_used|Mem_requested|Mem_used|CPUtime|CPUtime_mins|Walltime_req|Walltime_used|Walltime_mins|JobFS_req|JobFS_used|Efficiency|Service_units(CPU_hours)|Queue|NCPUS/task|
|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
|fastqc.o|	5|	5|	20.0GB|	20.0GB|	01:44:24|	104.40|	01:00:00|	00:21:30|	21.50|	100.0MB|	8.05MB|	0.97|	3.58|normal|1|
|bbduktrim.o|	5|	5|	80.0GB|	74.0GB|	04:46:22|	286.37|	03:00:00|	01:03:02|	63.03|	100.0MB|	8.05MB|	0.91|	42.02|normal|1|


# Supplementary scripts

* `sh multiqc_single_dataset.sh Batch_1` performs multiqc for a single sequencing batch, e.g. Batch_1, as specified by DATASET in ..<cohort>.config. Run after fastQC.

# Cite us to support us!
 
The RNASeq-DE pipeline can be cited as DOI: https://doi.org/10.48546/workflowhub.workflow.152.1  
  
If you use our pipelines, please cite us:  
  
Sydney Informatics Hub, Core Research Facilities, University of Sydney, 2021, The Sydney Informatics Hub Bioinformatics Repository, <date accessed>, https://github.com/Sydney-Informatics-Hub/Bioinformatics
