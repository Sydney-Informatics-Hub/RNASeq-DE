# RNASeq-DE

  - [Description](#description)
  - [User guide](#user-guide)
      - [Set up](#set-up)
      - [Running the pipeline](#running-the-pipeline)  
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
  - [References](#References)
  
# Description

RNASeq-DE is a __highly scalable__ workflow that pre-processes Illumina RNA sequencing data for differential expression (raw FASTQ to counts) on the __National Compute Infrastructure, Gadi__. The workflow was designed to efficiently process and manage large scale projects (100s of samples sequenced over multiple batches). 

In summary, the steps of this workflow include:

0. __Set up__
1. __QC of raw FASTQs__: FastQC and MultiQC to obtain quality reports on raw fastq files
2. __Trim raw FASTQs__: BBduk trim to trim 3' adapters and poly A tails. [Optional] - QC trimmed FASTQs.
3. __Mapping__: STAR for spliced aware alignment of RNA sequencing reads (FASTQ) to a reference genome
    * Prepare reference: STAR indexing for a given reference genome and corresponding annotation file
    * Perform mapping with trimmed FASTQs to prepared reference
    * Compress unmapped reads with pigz (optional)
    * Outputs: BAM per FASTQ pair (sequencing batch), unmapped reads, alignment stats
4. __Merge lane level to sample level BAMs__: SAMtools to merge and index sample BAMs
    * Merge multiplexed BAMs into sample level BAMs. For samples which are not multiplexed, files are renamed/symlinked for consistent file naming.
    * Outputs: `<sampleid>.final.bam` and `<sampleid>_<flowcell>.final.bam` (and index files)
5. __Mapping metrics__
    * RSeQC infer_experiment.py - check strand awareness, required for counting. Output: per sample reports which are summarized by cohort in `../QC_reports/<cohort>_final_bams_inter_experiment/<cohort>.txt`
    * RSeQC read_distribution.py - checks which features reads aligned to for each sample (summarized with multiqc). Expect ~30% of reads to align to CDS exons (provides total reads, total tags, total tags assigned. Groups by: CDS exons, 5' UTR, 3' UTR, Introns, TSS up and down regions). 
    * [Optional] RSeQC bam_stat.py - for each BAM, print QC metrics (numbers are read counts) including: Total records, QC failed, PCR dups, Non primary hits, unmapped, mapq, etc (similar metrics are provided by STAR, but can be used on sample level BAMs).
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

The scripts in this repository use relative paths and require careful setup to ensure that scripts can locate input files seamlessly. To start:

```
git clone https://github.com/Sydney-Informatics-Hub/RNASeq-DE
cd RNASeq-DE
```

### Required inputs and directory structure

Please provide the following files to run the workflow:

* __Illumina raw FASTQ files__ 
    * Organised into dataset/sequencing batch directories. Most sequencing companies will provide FASTQs in this structure. No other directory structure is supported.
    * FASTQ sequence identifier must follow the [standard Illumina format](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm) `@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>`. 
* __Reference files__
    *  Reference genome primary assembly (.fasta) and corresponding annotation (.gtf) file needs to be in a sub-directory in `Reference`. 
    *  Annotation in BED format is required for RSeQC's infer_experiment.py (CLI tools such as [gtf2bed.pl](https://github.com/ExpressionAnalysis/ea-utils/blob/master/clipper/gtf2bed) can convert GTF to BED format). 
    *  References can be obtained:
      * following recommendations in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). 
      * from [Ensembl](https://asia.ensembl.org/info/data/ftp/index.html)
* __.config file__: created using the [guide](#cohortconfig) below.

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

#### .config

The `.config` file is tab-delimited text file that is used to tell the scripts which samples to process, how to process them, and where it can locate relevant input files. An example is provided below:

|#FASTQ|	SAMPLEID|	DATASET|	REFERENCE|	SEQUENCING_CENTRE|	PLATFORM|	RUN_TYPE_SINGLE_PAIRED|	LIBRARY|
|------|---------|----------|------------------------|--------------------|-----------|-------------------------|--------|
|sample1_1.fastq.gz|	SAMPLEID1|	Batch_1|	GRCh38|	KCCG|	ILLUMINA|	PAIRED|	1|
|sample1_2.fastq.gz|     SAMPLEID1|       Batch_1| GRCh38|  KCCG|    ILLUMINA|        PAIRED|  1|
|sample2.fastq.gz|	SAMPLEID2|	Batch_2|	GRCm38|	KCCG|	ILLUMINA|	SINGLE|	1|

To create a `.config` using excel:

* Use column descriptions below to help you populate your config file
    * __All columns are required in this order and format__
* Save your `.config` file in the `RNASeq-DE` directory on Gadi.
    * Either copy and paste the contents into a text editor in Gadi, or save as a tab-delimited text file in excel. 
    * The file must be suffixed with `.config` 
    * File prefix: This is used to name outputs. Use something meaningful (e.g. your cohort name) and avoid whitespace.
    * Header lines must start with `#` 
   
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

Run all scripts from the `Scripts` directory once you have completed [set up](#set-up).

Generally, steps involve:
1. Running a `<task>_make_input.sh` script to prepare for parallel processing.
   * This makes an inputs file, e.g. `./Inputs/<task>.inputs`
   * This will often use your `<cohort>.config` file to know which files or samples you would like to process in a single job
2. Running a `<task>_run_parallel.pbs` script. 
   * This launches multiple tasks (e.g. `./Scripts/<task>.sh`) in parallel
   * Each line of `./Inputs/<task>.inputs` is used as input into a single `<task>.sh`
   * __Compute resources__ should be scaled to the size of your data. Recommendations are provided in the user guide.

The steps below explain how to process samples in `cohort.config`. __Replace `cohort.config` with the path to your own `.config` file.__
   
### 1. QC of raw FASTQs

This step performs FastQC to obtain quality reports per input FASTQ file. For multiple FASTQ files, reports can be summarized with MultiQC.

* Required inputs: `cohort.config`
   * "DATASET" is used locate file in "FASTQ" and name output directories
* Outputs: FastQC and MultiQC reports are written to `../dataset_fastQC`
   
To run FastQC for all raw FASTQ files in your `cohort.config` file, create input file for parallel processing:

```
sh fastqc_make_input.sh cohort.config
```

Edit `fastqc_run_parallel.pbs` by:
   * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
   * Adjusting PBS directive compute requests, scaling to your input size
      * Each `fastqc.sh` task requires NCPUS=1, 4 GB mem and ~00:30:00 walltime to process one FASTQ file with ~90 M reads. Scale walltime to number of expected reads per FASTQ.
      * For ~100 FASTQ files (or 50 FASTQ pairs), I recommend `-l walltime=01:30:00,ncpus=48,mem=190GB,wd`, `-q normal`. This will allow processing of 48 tasks in parallel.

Submit `fastqc_run_parallel.pbs` to perform FastQC in parallel (1 fastq file = 1 `fastqc.sh` task) by:
   
```
qsub fastqc_run_parallel.pbs
```

Once `fastqc_run_parallel.pbs` is complete, you can summarize reports using:

```
sh multiqc.sh ../dataset_fastQC
```

### 2. Trim raw FASTQs

This step trims raw FASTQ files using BBDuk trim. 

* Required inputs: `cohort.config`
   * "DATASET" is used locate file in "FASTQ" and name output directories
   * "RUN_TYPE_SINGLE_PAIRED" is used to indicate whether to trim as single or paired reads
* Outputs: Directory `../<dataset>_trimmed` containing trimmed FASTQ files
   
Task scripts `bbduk_trim_paired.sh` and `bbduk_trim_single.sh`are apply the following settings by default:
   *  Recommendated parameters under the "Adapter Trimming" example on the [BBDuk Guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) are used
   * `trimpolya=${readlen}`, where `${readlen}` is the length of your sequencing reads, obtained from the FASTQ file (assumes all read lengths in a single FASTQ are equal)
   * NO quality trimming is applied
   
To run BBDuk trim for FASTQ files in `cohort.config`, create input file for parallel processing:

```
sh bbduk_trim_make_input.sh cohort.config
```

Edit `bbduk_trim_run_parallel.pbs` by:
   * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
   * Adjusting PBS directive compute requests, scaling to your input size 
      * Each `bbduk_trim_paired.sh` task requires NCPU=6, 23 GB mem and ~00:19:00 walltime for 90 M pairs of FASTQ reads.
      * For ~15 FASTQ pairs, 90 M pairs each, I recommend: `-l walltime=02:00:00,ncpus=48,mem=190GB,wd`, `-q normal`

Submit `bbduk_trim_run_parallel.pbs`. This launches parallel tasks for: 1 FASTQ pair = 1 input for `bbduk_trim_paired.sh`, 1 FASTQ file = 1 input for `bbduk_trim_single.sh`:
   
```
qsub bbduk_trim_run_parallel.pbs
```

#### QC of trimmed FASTQs (optional)

You can check the quality of the data after trimming using:
  
```
sh fastqc_trimmed_make_input.sh cohort.config
```

Edit `fastqc_run_parallel.pbs` by:
   * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
   * Adjusting PBS directive compute requests, scaling to your input size (use your previous fastqc job as a guide)

```
qsub fastqc_trimmed_run_parallel.pbs
```

### 3. Mapping

#### Preparing your reference for STAR

Each reference under the "REFERENCE" column in `cohort.config` needs to be prepared with STAR before mapping. For most, this will only be one reference genome, prepared once per project. 

* Required inputs: reference genome in FASTA format (e.g. in `./Reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa`) and annotation file in GTF format (e.g. `./Reference/GRCh38/Homo_sapiens.GRCh38.103.gtf`). The "REFERENCE" column in your cohort.config file is used to locate the subdirectory (GRCh38 in this example), so make sure they match!
* Required parameters: overhang (read length - 1). The default is 149 for 150 base pair reads.

Edit the variables with the required inputs and parameters in `star_index.pbs` in the section shown below:
  
```
# Change dir, ref, gtf and overhang variables below
dir=../Reference/GRCh38
# ref and gtf files under dir
ref=Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=Homo_sapiens.GRCh38.103.gtf
# sjdbOverhang = ReadLength-1 (for creating splice junction database)
overhang=149
```
The default compute parameters are sufficient for human, mouse or other similar genome.

Submit the job:

```
qsub star_index.pbs
```

#### Mapping trimmed reads to prepared reference

This step will map trimmed FASTQ files to a prepared reference genome using STAR.

* Required inputs: `cohort.config`, `../<dataset>_trimmed` containing trimmed FASTQ files
    * "SAMPLEID" is used locate trimmed FASTQ in `../<dataset>_trimmed` and name output files
    * "PLATFORM", "SEQUENCING_CENTRE" from `.config`, and flowcell and lane from FASTQ sequence identifiers are used in BAM read group headers. 
    * "RUN_TYPE_SINGLE_PAIRED" is used to indicate whether to map as single or paired reads
* Outputs: Directory and output files prefixed `../<dataset>_STAR/${sampleid}_${lane}_`

To map all trimmed reads for to references specified in `cohort.config` file, prepare inputs for parallel processing by:
  
```
sh star_align_trimmed_make_input.sh cohort.config
```
  
`star_align_run_parallel.pbs` run task scripts `star_align_paired.sh` and/or `star_align_single.sh` and by default:
  * Will output coordinate sorted BAMs
  * Will output unmapped reads in FASTQ format (as pairs, if run type was paired)
  
Edit `star_align_run_parallel.pbs` by:
   * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
   * Adjusting PBS directive compute requests, scaling to your input size
        * Each `star_align_paired_with_unmapped.sh` task requires NCPUS=24, 96 GB mem and ~00:10:00 walltime to map 90 M pairs of FASTQ reads.
        * For ~15 FASTQ pairs, 90 M pairs each, I recommend: `-l walltime=01:30:00,ncpus=240,mem=950GB,wd`, `-q normal`
  
Submit `star_align_run_parallel.pbs`. This launches parallel tasks for: 1 trimmed FASTQ pair = 1 input for `star_align_paired_with_unmapped.sh`, 1 trimmed FASTQ file = 1 input for `star_align_single_with_unmapped.sh`: 

```
qsub star_align_run_parallel.pbs
```
  
#### Compress unmapped reads with pigz (optional)
  
STAR outputs unmapped reads as unzipped FASTQ files. To save disk, we can compress these files using pigz. 

* Required inputs: STAR unmapped FASTQ files. Filenames end in "*Unmapped.out.mate1", "*Unmapped.out.mate2"
* Outputs: STAR unmapped gzipped FASTQ files. Filenames end in "*Unmapped.out.mate1.gz", "*Unmapped.out.mate2.gz"
  
To do this, edit `pigz_run_parallel.pbs` by:
  * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
  * This script also creates inputs for parallel processing. By default, it will search for all unmapped STAR files using STAR's file naming convention. You may want to check the find commands on the command line for your version of STAR.
  *  Adjusting PBS directive compute requests, scaling to your input size
    * For ~200 unmapped files (or ~100 pairs), I suggest `-l walltime=02:00:00,ncpus=28,mem=128GB,wd`, `-q normalbw`

The task `pigz.sh` will delete the original file by default.
  
Submit your job by:

```
qsub pigz_run_parallel.pbs
```
  
### 4. Merging lane level BAMs into sample level BAMs

This step merges sample lane level BAMs into sample level BAMs (skipped if samples were not multiplexed). All final BAMs are organised into `cohort_final_bam` directory and are then indexed with SAMtools.

* Required inputs: `cohort.config` and sample BAMs in `*STAR` directories
  * A unique list of sample IDs are taken from `cohort.config`
* Outputs: Final BAMs in `cohort_final_bam/<sampleid>.final.bam` and index files `cohort_final_bam/<sampleid>.final.bam.bai`
  * For samples that were not multiplexed, STAR bams are symlinked into the `cohort_final_bam` directory.
  * For multiplexed samples, a sample level and flowcell level (symlinked) BAM will be created in the `cohort_final_bam` directory. Flowcell level BAMs are useful for checking technical batch effects.

To obtain final BAMs for sample IDs in `cohort.config`, create input file for parallel processing:

```
sh samtools_merge_index_make_input.sh cohort.config
```

Edit `samtools_merge_index_run_parallel.pbs` by:
  * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
  * Adjusting PBS directive compute requests, scaling to your input size 
      * This will be highly dependant on how many samples in your cohort are multiplexed
      * For multiplexed samples, allow NCPUS=3
      * For non-multiplexed samples, allow NCPUS=1
      * For ~882 samples with an average of 80 M paired reads, including ~10 multiplexed, I suggest `-l walltime=5:00:00,ncpus=48,mem=190GB,wd`, `-q normal`

Submit the job by:

```
qsub samtools_merge_index_run_parallel.pbs
```

### 5. Mapping metrics

#### RSeQC's infer_experiment.py

This step uses RSeQC's infer_experiment.py to infer the library strand awareness (i.e. forward, reverse, or not strand aware) that was used to prepare the samples for sequencing. This is required for `htseq-count`. 
  
* Required inputs: A directory of BAM files (e.g. `cohort_final_bams`) and a reference annotation file in BED format (CLI tools such as [gtf2bed](https://github.com/ExpressionAnalysis/ea-utils/blob/master/clipper/gtf2bed) can convert GTF to BED format.
* Outputs: 
    * infer_experiment.py results per BAM as <sampleid>.txt or <sampleid>_<flowcell>.txt files
    * One summary table for all BAMs for paired data in ../QC_reports/${outfileprefix}_infer_experiment/cohort_final_bams.txt with `#FILE REVERSE FORWARD` headers. REVERSE is the proportion of reads supporting reverse strand awareness (1-+,2++,2--) and FORWARD is the proportion of reads supporting forward strandawareness (1++,1--,2+-,2-+). 
  
The `infer_experiment_final_bams.sh` script processes multiple BAMs in parallel on the login node/command line. With the path to the directory containing `*final.bam` files, e.g. ../cohort_final_bams:
  
```
sh infer_experiment_final_bams.sh ../cohort_final_bams
```
  
#### RSeQC's read_distribution.py

This step uses RSeQC's read_distribution.py to check the distribution of aligned reads across features (e.g. exon, intron, UTR, intergenic, etc). 
  
  * Required inputs: `cohort.config`, <sampleid>.final.bam and reference annotation file in BED format
  * Output: Per sample output in `../QC_reports/<cohort>_read_distribution/<sampleid>_read_distribution.txt`
  
To obtain `read_distribution.py` reports for sample BAMs in `cohort.config`, create input file for parallel processing:

```
sh read_distribution_make_input.sh cohort.config
```
`read_distribution_run_parallel.pbs` will run task script `read_distribution.sh` with `read_distribution.py` default settings applied.
 
 Edit `read_distribution_run_parallel.pbs` by:
  * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
  * Adjusting PBS directive compute requests, scaling to your input size
      * Each task will only require NCPUS=1, ~00:25:00 walltime (1 BAM with ~80 M paired reads)
      * For ~260 BAM files, ~80 M pairs of reads per sample, I suggest: `-l walltime=04:00:00,ncpus=56,mem=256GB,wd`, `-q normalbw`
 
Once `read_distribution_run_parallel.pbs` is complete, you can summarize reports using:

```
sh multiqc.sh ../QC_reports/cohort_read_distribution
```
  
#### RSeQC's bam_stat.py (optional)
  
This step uses RSeQC's bam_stat.py to check alignment metrics of BAM files, including: Total records, QC failed, PCR dups, Non primary hits, unmapped, mapq, etc. 
  
  * Required inputs: `cohort.config` and `cohort_final_bams/<sampleid>.final.bam` files
  * Output: Per sample output in `../QC_reports/<cohort>_final_bams_bam_stat/<sampleid>_bam_statc.txt`

To obtain `bam_stat.py` reports for sample BAMs in `cohort.config`, create input file for parallel processing:

```
sh bam_stat_make_input.sh cohort.config
```

`bam_stat_run_parallel.pbs` will run task script `bam_stat.sh` with `bam_stat.py` default settings applied.

Edit `bam_stat_run_parallel.pbs` by:
  * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
  * Adjusting PBS directive compute requests, scaling to your input size
      * Each task requires NCPUS=1
      * For 907 BAMs with ~80 M pairs of reads, I suggest `-l walltime=04:00:00,ncpus=112,mem=512GB,wd`, `-q normalbw`

Once `bam_stat_run_parallel.pbs` is complete, you can summarize reports using:

```
sh multiqc.sh ../QC_reports/cohort_final_bams_bam_stat
```

#### summarise_STAR_alignment_stats.pl
 
The `summarise_STAR_alignment_stats.pl` script will collate mapping information from STAR output for each aligned file. 

* Required inputs: `dataset_STAR/<sampleid>_<lane>_Log.final.out` files, obtained from mapping
* Output: `QC_reports/cohort_STAR_metrics.txt`

Run the following command on the login node, on the command line:

```
perl summarise_STAR_alignment_stats.pl cohort.config
```
#### SAMtools idxstats

This step runs samtools idxstats for all BAMs in a directory, and summarize the reports with multiQC.

* Required inputs: Directory with `.bam` files
* Outputs: `QC_reports/<outfileprefix>_samtools_idxstats`

Run this on the login node, providing the path to your directory containing BAM files, e.g.:
  
```
sh samtools_idxstats_final_bams.sh ../cohort_final_bams
```

### 6. Raw counts

This step uses htseq-count to obtain aligned read counts present across features in a genome.

#### Counts per sample
 
* Required inputs: `cohort.config`, `cohort_final_bams/<sampleid>.final.bam` and `Reference/GRCh38/*gtf` file
* Output: Per sample counts in `cohort_htseq-count/<sampleid>.counts`
  
To obtain counts from final BAMs for sample IDs in `cohort.config`:
  * Edit the strand= variable in `htseq-count_custom_make_input.sh` if your libraries are not reverse strand aware!
  * Create input file for parallel processing:

```
sh htseq-count_custom_make_input.sh cohort.config
```
Note: `htseq-count_custom_make_input.sh` will search for all .final.bam in cohort_final_bams, including sample flowcell level BAMs if available. To use only sample level bams, create inputs with `htseq-count_make_input.sh`

 `htseq-count_run_parallel.pbs` will run task script `htseq-count.sh` with htseq-count recommended settings.
  
Edit `htseq-count_run_parallel.pbs` by:
  * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
  * Adjusting PBS directive compute requests, scaling to your input size
      * Each task requires NCPUS=1. Walltime will scale to the number of reads in a BAM file. A sample with 130 M paired reads requires ~04:00:00 walltime. A sample with 80 M paired reads requires ~02:20:00 walltime.

#### Cohort count matrix (optional)

This step creates a standard count matrix file (genes = rows, sampleIDs = columns) using htseq-count output. All samples within your `.config` file that have `<sampleid>.counts` file available will be included in the final matrix. 
  
* Required inputs: `cohort.config`, used to locate `cohort_htseq-count` directory and <sampleid>.counts
* Output: Count matrix file in `<cohort>_htseq-count/<cohort>.counts`

For smaller cohorts (<100 samples), run on the login node:

```
perl htseq-count_make_matrix_custom.pl cohort.config
```

For very large cohorts (>100 samples), run as a job by editing `htseq-count_make_matrix_custom.pbs` by:
  * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
  * Providing the `config=` variable the path to your cohort.config file
  * For most cohorts, the default compute resources should be sufficient. Larger cohorts will require more memory.

Submit the job by:

```
qsub htseq-count_make_matrix_custom.pbs
```

### 7. Normalize counts

This step quantifies transcript abundance for genomic features (gene, transcript, exon, intron) as transcripts per million (TPM) normalized counts. 

#### TPM counts per sample
  
* Required inputs: `cohort.config`, `Reference/GRCh38/*gtf` and `cohort_final_bams` directory containing indexed BAMs
* Outputs: TPMCalculator outputs per input BAM in directory `cohort_TPMCalculator`
  
To obtain TPM normalized counts across features for samples in `cohort.config`, create input file by:
* Providing the `gtf=` variable the path to your GTF file (by default this is `../Reference/GRCh38/Homo_sapiens.GRCh38.103.gtf`
* Create the input file for parallel processing:

```
sh tpmcalculator_make_input.sh cohort.config
```

`tpmcalculator_run_parallel.pbs` will run task script `tpmcalculator.sh`. By default, this applies the TPMCalculator's developer's settings and addtional settings including:
  * `-a` to print all features. Required for the next step, which collates sample TPM counts into a cohort count matrix
  * `-p` use only properly paired reads. Recommended by the TPMCalculator developer
  * `-e` extended output, to include transcript level TPM values
  * `-q 255` apply minimum MAPQ value of 255 to filter reads. This value is [recommended by the developer](https://github.com/ncbi/TPMCalculator/issues/72) for STAR aligned BAMs 
 
Edit `tpmcalculator_run_parallel.pbs` by:
  * Replacing PBS directive parameters, specifically <project> with your NCI Gadi project short code
  * Adjusting PBS directive compute requests, scaling to your input size
    * Each task requires NCPUS=1, 3-4 GB memory
    * For 907 BAMs with ~80 M paired reads each, I suggest `-l walltime=8:00:00,ncpus=768,mem=3040GB,wd`, `-q normal`
 
#### TPM cohort count matrix (optional)

This step creates two TPM cohort count matricies, one at the gene level, the other at the transcript level.
  
  * Required inputs: Directory containing sample level TPMCalculator outputs
  * Output: `TPM_GeneLevel.counts` and `TPM_TranscriptLEvel.counts` in the input directory provided

Edit the `tpmcalculator_make_matrix.pbs` script:
  * Providing the `tpmdir=` variable the path to your TPMCalculator directory containing TPMCalculator outputs per sample
  * Adjust memory if required. 54 GB memory was required for ~907 samples

Run the script:

```
qsub tpmcalculator_make_matrix.pbs
```
         
# Benchmarking

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
|Version            | 1.0.0                 |
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

The software listed below are used in the RNASeq-DE pipeline. Some of these are installed globally on NCI Gadi (check with `module avail` for the current software). Install python3 packages by `module load python3/3.8.5`, and then using the `pip3 install` commands. These will be installed in `$HOME`. All other software need to be installed in your project's `/scratch` directory and module loadable. 

openmpi/4.0.2 (installed globally)

nci-parallel/1.0.0a (installed globally)

SAMtools/1.10 (installed globally)

python3/3.8.5 (installed globally)

fastQC/0.11.7

multiqc/1.9

BBDuk/37.98

STAR/2.7.3a

rseqc/4.0.0

htseq-count/0.12.4

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

## Authors 

- Tracy Chew (Sydney Informatics Hub, University of Sydney)
- Rosemarie Sadsad

## Acknowledgements 
Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities. We suggest including the following acknowledgement in any publications that follow from this work:

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia.

Documentation was created following the [Australian BioCommons documentation guidelines](https://github.com/AustralianBioCommons/doc_guidelines).

## Cite us to support us!

 Chew, T., & Sadsad, R. (2022). RNASeq-DE (Version 1.0) [Computer software]. https://doi.org/10.48546/workflowhub.workflow.152.1
 
# References

Anders, S., Pyl, P.T., Huber, W., 2014. HTSeq--a Python framework to work with high-throughput sequencing data. Bioinformatics. https://doi.org/10.1093/bioinformatics/btu638

Andrews, S. 2010. FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]
  
Aronesty, E. 2011. ea-utils : "Command-line tools for processing biological sequencing data"; https://github.com/ExpressionAnalysis/ea-utils

BBMap - Bushnell B. - sourceforge.net/projects/bbmap/

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., Gingeras, T.R., 2012. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. https://doi.org/10.1093/bioinformatics/bts635
  
Ewels, P., Magnusson, M., Lundin, S., Käller, M., 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. https://doi.org/10.1093/bioinformatics/btw354
  
Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R., & 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics (Oxford, England), 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

Vera Alvarez, R., Pongor, L.S., Mariño-Ramírez, L., Landsman, D., 2018. TPMCalculator: one-step software to quantify mRNA abundance of genomic features. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty896
  
Wang, L., Wang, S., Li, W., 2012. RSeQC: quality control of RNA-seq experiments. Bioinformatics. https://doi.org/10.1093/bioinformatics/bts356
