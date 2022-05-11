# RNASeq-DE

  - [Description](#description)
  - [User guide](#user-guide)
  - [Benchmarking](#benchmarking)
  - [Workflow summaries](#workflow-summaries)
      - [Metadata](#metadata)
      - [Component tools](#component-tools)
      - [Required (minimum)
        inputs/parameters](#required-minimum-inputsparameters)
  - [Help/FAQ/Troubleshooting](#helpfaqtroubleshooting)
  - [3rd party Tutorials](#3rd-party-tutorials)
  - [Licence(s)](#licences)
  - [Acknowledgements/citations/credits](#acknowledgementscitationscredits)

# Description

RNASeq-DE is a highly scalable workflow that pre-processes RNA sequencing data for differential expression (raw FASTQ to counts) on the __National Compute Infrastructure, Gadi__. The workflow was designed to efficiently process and manage large scale projects (100s of samples sequenced over multiple batches). 

In summary, the steps of this workflow include:

0. __Set up__
1. __QC of raw FASTQs__: FastQC and MultiQC to obtain quality reports on raw fastq files
2. __Trim raw FASTQs__: BBduk trim to trim 3' adapters and poly A tails. [Optional] - QC trimmed FASTQs.
3. __Mapping__: STAR for spliced aware alignment of RNA sequencing reads (FASTQ) to a reference genome
    * Prepare reference: STAR indexing for a given reference genome and corresponding annotation file
    * Perform mapping with trimmed FASTQs to prepared reference
    * Pigz is used to gzip unmapped files
    * Outputs: BAM per FASTQ pair (sequencing batch), unmapped reads, alignment stats
4. __Merge lane level to sample level BAMs__: SAMtools to merge and index sample BAMs
    * Merge multiplexed BAMs into sample level BAMs. For samples which are not multiplexed, files are renamed/symlinked for consistent file naming.
    * Outputs: `<sampleid>.final.bam` and `<sampleid>_<flowcell>.final.bam` (and index files)
5. __QC of mapping__
    * RSeQC infer_experiment.py - check strand awareness, required for counting. Output: per sample reports which are summarized by cohort in `../QC_reports/<cohort>_final_bams_inter_experiment/<cohort>.txt`
    * [Optional] RSeQC bam_stat.py - for each BAM, print QC metrics (numbers are read counts) including: Total records, QC failed, PCR dups, Non primary hits, unmapped, mapq, etc (similar metrics are provided by STAR, but can be used on sample level BAMs).
    * RSeQC read_distribution.py - checks which features reads aligned to for each sample (summarized with multiqc). Expect ~30% of reads to align to CDS exons (provides total reads, total tags, total tags assigned. Groups by: CDS exons, 5' UTR, 3' UTR, Introns, TSS up and down regions). 
    * [Optional] `summarize_STAR_alignment_stats.pl`: collates per sample STAR metric per flowcell level BAM (use read_distribution for sample level BAMs). Uses datasets present in a cohort.config file to find these BAMs. Inputs: per sample `*Log.final.out`. Output: `../QC_reports/<cohort>_STAR_metrics.txt
    * [Optional] SAMtools idxstats: summarize number of reads per chromosome (useful for inferring gender, probably needs a bit more work)
6. __Raw counts__: HTSeq
    * Count reads in `<sampleid>.final.bam` that align to features in your reference and create a count matrix for a cohort
    * Output: `<sample>.counts` per input BAM and a count matrix as `<cohort>.counts`
7. __Normalized counts__: BAMtools/TPMCalculator
    * Obtain TPM normalized counts for gene and transcripts in your reference from `<sampleid>.final.bam` and create a TPM count matrix for a cohort
    * Output: Per sample TPM counts and cohort count matricies as `TPM_TranscriptLevel.counts`, `TPM_GeneLevel.counts`

# User guide

## Set up

The scripts in this repository use relative paths and require careful setup to ensure that scripts can locate input files seamlessly. To start, go to your working directory (e.g. /scratch/ab1 for your ab1 project on NCI Gadi) and:

```
git clone https://github.com/Sydney-Informatics-Hub/RNASeq-DE
cd RNASeq-DE
```

### Required inputs and directory structure

Please provide the following files to run the workflow:

* __Raw FASTQ files__ organised into dataset directories. Most sequencing companies will provide FASTQs in this structure. 
* __Reference files__. Reference genome primary assembly (.fasta) and corresponding annotation (.gtf) file needs to be in a sub-directory in `Reference`. References can be obtained:
    * following recommendations in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). 
    * from [Ensembl](https://asia.ensembl.org/info/data/ftp/index.html)
* __.config file__: create using the [guide](#cohortconfig) below.

Your __RNASeq-DE__ directory structure should match the following:

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
├── cohort.config
└── Scripts
```

### cohort.config

The `cohort.config` file is used to tell the scripts which samples to process, how to process them, and where it can locate relevant input files. An example is provided below:

|#FASTQ|	SAMPLEID|	DATASET|	REFERENCE|	SEQUENCING_CENTRE|	PLATFORM|	RUN_TYPE_SINGLE_PAIRED|	LIBRARY|
|------|---------|----------|------------------------|--------------------|-----------|-------------------------|--------|
|sample1_1.fastq.gz|	SAMPLEID1|	Batch_1|	GRCh38|	KCCG|	ILLUMINA|	PAIRED|	1|
|sample1_2.fastq.gz|     SAMPLEID1|       Batch_1| GRCh38|  KCCG|    ILLUMINA|        PAIRED|  1|
|sample2.fastq.gz|	SAMPLEID2|	Batch_2|	GRCm38|	KCCG|	ILLUMINA|	SINGLE|	1|

To create a `cohort.config file` using excel:

* Open or copy template headers into excel. `cohort.config` is a tab-delimited text file
* Use column descriptions below to help you populate your config file
   * __All columns are required in this order and format__
* Save your `cohort.config` file in the `RNASeq-DE` directory. 
   * `cohort` is used to name output files and directories. Change the prefix of this file to something more meaningful.
   
Column descriptions for __cohort.config__:

|Column name| Description|
|----|--------|
|#FASTQ| FASTQ filename. This column can be populated with `ls *f*q.gz | sort -V` in your sequencing batch directory. Paired files are recognised by conventional filenames - 1 or 2 at the end of the filename, before the fastq.gz extension e.g. `<sampleid>_<flowcell>_<tag>_<lane>_R<1|2>.fastq.gz`|
|SAMPLEID| The sample identifier used in your laboratory. Sample IDs are used to name output files. Avoid whitespace.|
|DATASET| Directory name containing the sample FASTQs. This is typically analogous to the sequencing batch that the FASTQ file was generated. |
|REFERENCE| Reference subdirectory name, e.g. GRCh38 or GRCm38 in the above example (must case match). Scripts will use reference files (.fasta and .gtf) and STAR index files for the FASTQ file/sample for alignment and counting. |
|SEQUENCING_CENTRE| e.g. AGRF. This is used in the read group header in the output BAM file for the aligned FASTQ. Avoid whitespace.|
|PLATFORM| e.g. ILLUMINA. This is used in the read group header in the output BAM file for the aligned FASTQ.|
|RUN_TYPE_SINGLE_PAIRED| Input SINGLE or PAIRED. This is used to indicate whether you want to process single read data or paired end data (STAR and BBduk trim).| 
|LIBRARY| The sequencing library of the FASTQ file. This is used in the read group header in the output BAM file for the aligned FASTQ. Use 1 if unknown. No whitespace please.|

## Running the pipeline

When you have completed [Set up](#set-up), change into the Scripts directory using `cd Scripts` and run all scripts from here. 

Generally, steps involve:
1. Running a `<task>_make_input.sh` script to prepare for parallel processing.
  * This makes an inputs file, e.g. `./Inputs/<task>.inputs`
  * This will often use your `<cohort>.config` file to know which files or samples you would like to process in parallel
2. Running a `<task>_run_parallel.pbs` script. 
  * This launches multiple tasks (e.g. `./Scripts/<task>.sh`) in parallel
  * Each line of `./Inputs/<task>.inputs` is used as input into a single `<task>.sh`
  * __Compute resources__ should be scaled to the size of your data, using this user guide or [benchmarking](#benchmarking) as a guide
   
### 1. QC of raw FASTQs

This step performs FastQC to obtain quality reports per input FASTQ file. For multiple FASTQ files, reports can be summarized with MultiQC.

* Required inputs: `cohort.config`
   * "DATASET" is used locate file in "FASTQ" and name output directories
* Outputs: FastQC and MultiQC reports are written to `../dataset_fastQC`
   
To run FastQC for all raw FASTQ files in your `cohort.config` file, create input file for parallel processing:

```
sh fastqc_make_input.sh cohort.config
```

Edit `qsub fastqc_run_parallel.pbs` by:
   * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
   * Adjusting PBS directive compute requests, scaling to your input size (consider number of tasks, size of FASTQ)
      * Each task requires NCPU=1. 1 FASTQ file containing ~15 M reads takes about 00:05:00 in walltime. 

Submit `qsub fastqc_run_parallel.pbs` to perform FastQC in parallel (1 fastq file = 1 `fastqc.sh` task) by:
   
```
qsub fastqc_run_parallel.pbs
```

Once `fastqc_run_parallel.pbs` is complete, you can summarize reports using:

```
sh multiqc.sh ../dataset_fastQC`
```

### 2. Trim raw FASTQs

This step trims raw FASTQ files using BBDuk trim. 

* Required inputs: `cohort.config`
   * "DATASET" is used locate file in "FASTQ" and name output directories
   * "RUN_TYPE_SINGLE_PAIRED" is used to indicate whether to trim as single or paired reads
* Outputs: Directory `../<dataset>_trimmed` containing trimmed FASTQ files
   
Task scripts `bbduk_trim_paired.sh` and `bbduk_trim_single.sh`are apply the following settings by default:
   *  Recommendated parameters under the "Adapter Trimming" example on the [BBDuk Guide](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) are used
   * `trimpolya=${readlen}`, where `${readlen}` is the length of your sequencing reads, obtained from the FASTQ file (assumes all read lengths in a single FASTQ are equal)
   * NO quality trimming is applied
   
To run BBDuk trim for all raw FASTQ files in your `cohort.config` file, create input file for parallel processing:

```
sh bbduk_trim_make_input.sh cohort.config
```

Edit `bbduk_trim_run_parallel.pbs` by:
   * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
   * Adjusting PBS directive compute requests, scaling to your input size (consider number of tasks, size of FASTQ)
      * NCPU=6

Submit `qsub bbduk_trim_run_parallel.pbs` to run bbduk.sh in parallel (1 FASTQ pair = 1 input for `bbduk_trim_paired.sh`, 1 FASTQ file = 1 input for `bbduk_trim_single.sh`):
   
```
qsub bbduk_trim_run_parallel.pbs
```

QC of trimmed FASTQs

### 3. Mapping

### 4. Merging lane level BAMs into sample level BAMs

### 5. Collect BAM QC metrics

### 6. Raw counts

### 7. Normalize counts
         
# Benchmarking

The below metrics were obtained for a human cohort comprised of:
* 96 FASTQ pairs
* ~80 M paired reads (each read 150 bp length)

| #JobName                              | CPUs_requested | CPUs_used | Mem_requested | Mem_used | CPUtime    | CPUtime_mins | Walltime_req | Walltime_used | Walltime_mins | JobFS_req | JobFS_used | Efficiency | Service_units(CPU_hours) | Job_exit_status | Date       | Time     |
|---------------------------------------|----------------|-----------|---------------|----------|------------|--------------|--------------|---------------|---------------|-----------|------------|------------|--------------------------|-----------------|------------|----------|
| bam_stat.o                            | 112            | 112       | 512.0GB       | 269.54GB | 208:45:39  | 12525.65     | 10:00:00     | 2:05:41       | 125.68        | 400.0MB   | 8.26MB     | 0.89       | 293.26                   | 0               | 23/02/2022 | 11:16:58 |
| bbduk_trim.o                          | 48             | 48        | 190.0GB       | 169.67GB | 22:36:20   | 1356.33      | 4:00:00      | 0:40:14       | 40.23         | 100.0MB   | 8.17MB     | 0.7        | 64.37                    | 0               | 22/02/2022 | 15:29:34 |
| fastqc.o                              | 30             | 30        | 150.0GB       | 106.41GB | 8:18:12    | 498.2        | 5:00:00      | 0:26:32       | 26.53         | 100.0MB   | 8.12MB     | 0.63       | 33.17                    | 0               | 22/02/2022 | 15:16:39 |
| fastqc_trimmed.o                      | 30             | 30        | 150.0GB       | 97.6GB   | 8:12:31    | 492.52       | 2:00:00      | 0:25:08       | 25.13         | 100.0MB   | 8.12MB     | 0.65       | 31.42                    | 0               | 22/02/2022 | 15:58:45 |
| htseq-count.o                         | 240            | 240       | 950.0GB       | 582.53GB | 2562:06:24 | 153726.4     | 30:00:00     | 19:51:09      | 1191.15       | 500.0MB   | 8.33MB     | 0.54       | 9529.2                   | 0               | 24/02/2022 | 4:58:30  |
| htseq-count_matrix.o                  | 1              | 1         | 32.0GB        | 16.02GB  | 0:02:53    | 2.88         | 5:00:00      | 0:03:45       | 3.75          | 100.0MB   | 0B         | 0.77       | 3                        | 0               | 24/02/2022 | 9:03:09  |
| multiqc_all_datasets_trimmed_fastQC.o | 1              | 1         | 32.0GB        | 26.59GB  | 0:06:07    | 6.12         | 5:00:00      | 0:09:37       | 9.62          | 100.0MB   | 1.45MB     | 0.64       | 7.69                     | 0               | 4/03/2022  | 15:28:28 |
| pigz.o                                | 28             | 28        | 128.0GB       | 33.98GB  | 1:25:20    | 85.33        | 2:00:00      | 0:07:15       | 7.25          | 100.0MB   | 8.1MB      | 0.42       | 4.23                     | 0               | 22/02/2022 | 17:32:47 |
| read_distribution.o                   | 144            | 144       | 570.0GB       | 419.27GB | 278:19:37  | 16699.62     | 10:00:00     | 2:28:54       | 148.9         | 300.0MB   | 8.4MB      | 0.78       | 714.72                   | 0               | 23/02/2022 | 11:42:52 |
| samtools_merge_index.o                | 48             | 48        | 190.0GB       | 99.15GB  | 86:30:00   | 5190         | 10:00:00     | 2:44:15       | 164.25        | 100.0MB   | 8.25MB     | 0.66       | 262.8                    | 0               | 22/02/2022 | 20:18:34 |
| star_align_trimmed_unmapped_out.o     | 240            | 240       | 950.0GB       | 736.62GB | 42:09:54   | 2529.9       | 24:00:00     | 0:21:28       | 21.47         | 500.0MB   | 8.17MB     | 0.49       | 171.73                   | 0               | 22/02/2022 | 15:59:35 |
| tpmcalculator_transcript.o            | 768            | 768       | 2.97TB        | 2.63TB   | 978:52:39  | 58732.65     | 10:00:00     | 5:14:53       | 314.88        | 1.56GB    | 8.38MB     | 0.24       | 8061.01                  | 0               | 23/02/2022 | 14:31:26 |
| tpmtranscript_matrix.o                | 1              | 1         | 96.0GB        | 53.19GB  | 0:37:32    | 37.53        | 5:00:00      | 0:45:33       | 45.55         | 100.0MB   | 0B         | 0.82       | 109.32                   | 0               | 24/02/2022 | 13:54:22 |

# Workflow summaries

## Metadata

|metadata field     | workflow_name / workflow_version  |
|-------------------|:---------------------------------:|
|Version            | 1.0                 |
|Maturity           | stable                            |
|Creators           | Tracy Chew                 |
|Source             | NA                                |
|License            | GNU GENERAL PUBLIC LICENSE        |
|Workflow manager   | None                          |
|Container          | None                              |
|Install method     | Manual                            |
|GitHub             | NA                                |
|bio.tools 	        | NA                                |
|BioContainers      | NA                                | 
|bioconda           | NA                                |


## Component tools

The software listed below are used in the RNASeq-DE pipeline. Some of these are installed globally on NCI Gadi (check with `module avail` for the current software). Install python3 packages by `module load python3/3.8.5`, and then using the `pip3 install` commands as listed below. These will be installed in `$HOME`. All other software need to be installed in your project's `/scratch` directory and module loadable. 

openmpi/4.0.2 (installed globally)

nci-parallel/1.0.0a (installed globally)

SAMtools/1.10 (installed globally)

python3/3.8.5 (installed globally)

pip3 install multiqc

pip3 install RSeQC

pip3 install HTSeq

fastQC/0.11.7

BBDuk/37.98

STAR/2.7.3a

rseqc/4.0.0

samtools/1.10

bamtools/2.5.1, TPMCalculator/0.0.4

## Required (minimum) inputs/parameters

* Short read FASTQ files (single or paired)
* Reference genome (FASTA), annotation (GTF) files. These can be obtained from [Ensembl FTP](http://www.ensembl.org/info/data/ftp/index.html)

# Help / FAQ / Troubleshooting

* Contact NCI for NCI related queries
* Contact tool developers for tool specific queries
* For RNASeq-DE workflow queries, please submit a [Github issue](https://github.com/Sydney-Informatics-Hub/RNASeq-DE/issues)

# License(s)

GNU General Public License v3.0

# Acknowledgements/citations/credits

### Authors 

- Tracy Chew (Sydney Informatics Hub, University of Sydney)
- Rosemarie Sadsad

### Acknowledgements 
Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities. We suggest including the following acknowledgement in any publications that follow from this work:

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia.

Documentation was created following the [Australian BioCommons documentation guidelines](https://github.com/AustralianBioCommons/doc_guidelines).

# Cite us to support us!

 Chew, T., & Sadsad, R. (2022). RNASeq-DE (Version 1.0) [Computer software]. https://doi.org/10.48546/workflowhub.workflow.152.1
