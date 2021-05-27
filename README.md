# Ultimate ATAC-seq Data Processing & Analysis Pipeline
From r**A**w (unaligned) BAM files to normali**Z**ed counts

A Snakemake implementation of the [BSF's](https://www.biomedical-sequencing.org/) [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline") extended by downstream processing and unsupervised analyses steps using bash, python and R. Reproducibility is ensured by using conda and Singularity(soon).

# Features
- Processing
    - alignment of both single-end and paired-end reads in raw/unaligned BAM format (bowtie2)
    - peak calling (macs2)
    - peak annotation and motif analysis (homer)
    - MultiQC report generation (multiqc)
- Quantification
    - consensus region set creation
    - consensus region set annotation (uropa with regulatory build and gencode; homer)
    - read count and peak support quantification of the consensus region set across samples
- optional: split data in multiple data sets (eg by cell type)
- Downstream Analysis (performed on the whole data set and each split separately)
    - region filtering
    - normalization by two different methods separately (TMM & CQN)
    - highly variable region (HVR) selection for each normalized data set
    - dimensionality reduction by PCA & UMAP and plotting of provided metadata after each step
    - mean-variance-relationship plotafter each step
- Snakemake report generation for workflow statistics, documentation, reproducibility and result inspection

# Installation (<10 minutes)
1. install snakemake, which requires conda & mamba, according to the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. clone/download this repository

# Configuration
You need 2 configuration files and 2 annotation files to run the complete workflow from A to Z:
- pipeline configuration (conifg/pipeline_config): one off set up in your environment
    - before execution edit the path to your project configuration (project_config)
    - if you are a [CeMMie](https://cemm.at/) you can just take the provided file
- project configuration: different for every project/dataset
- sample annotation: technical units
- metadata annotation: for downstream analysis steps

detailed descriptions can be found in the appendix below

# Execution
## Change working directory & activate conda environment
Execution always from within top level of the workflow directory (ie atacseq_pipeline/).
Snakemake commands only work from within the snakemake conda environment.
```
cd atacseq_pipeline
conda activate snakemake
```
## Execute a dry-run
command for a dry-run with option -n (-p makes Snakemake print the resulting shell command for illustration)
```
snakemake -p -n
```
## Execute workflow local or on a cluster
### Local execution
command for execution with one core
```
snakemake -p -j1 --use-conda
```
### Cluster execution
command for **vanilla cluster execution** on cluster engines that support shell scripts and have access to a common filesystem, (e.g. the Sun Grid Engine), more info in the [Snakemake Cluster Execution documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
```
snakemake -p --use-conda --cluster qsub -j 32
```

command for **cluster execution by using --profile**, submits every task as separate job with dependencies
```
snakemake -p --use-conda --profile config/slurm.cemm
```
the profile for CeMM's slurm environment is provided in the config/ directory, of note: 
- the number of jobs in the slurm.cemm/config.yaml should be set as high as necessary, because arrayed job subsmission does not work (yet) and the scheduler (eg SLURM) should take care of the priorization
- jobs which dependencies can never be fulfilled are automatically removed from the queue

If you are using another setup get your cluster execution profile here: [The Snakemake-Profiles project](https://github.com/snakemake-profiles/doc)

## Singularity execution (soon)
command for execution with singularity with flag --use-singularity
```
snakemake -p --use-conda --use-singularity
```
# Use as module in another Snakemake workflow (soon)
- [https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules)
- [https://slides.com/johanneskoester/snakemake-6#/8](https://slides.com/johanneskoester/snakemake-6#/8)

# Report
The pipeline automatically generates a self contained HTML based report in a zip archive containing the following sections:
- Workflow: interactive rulegraph to recapitulate individual steps, used software and conrete code (reproducibility)
- Statistics: duration and timing of individual steps
- Configuration: used pipeline configuration (accountability)
- Results
    - Configuration: all 3 used configuration files (project, samples, metadata)
    - QC reports: link to the MultiQC report
    - all: each step of the downstream analysis has its own searchable section with step-specific and unsupervised analysis plots
    - optional: one additional section per split, containing the respective downstream analysis results
    
# Results
Project directory structure:
- all: Downstream analysis results of the whole data, including consensus region set and region annotation
- atacseq_hub: genome browser track files (.bigWig) for each sample
- atacseq_report: statistics and metrics from the processing part as input for the MultiQC report
- atacseq_results: one directory per sample with all the processing and quantification results
- projectName_report.zip: self contained HTML Snakemake report
- split1: Downstream analysis results of subset 1 of the whole data
- split2: Downstream analysis results of subset 2 of the whole data

# Examples
We provide configuration files for two example datasets (mm10 & hg38).
Additionally, the report zip archive of the hg38 test example is provided to showcase the pipeline results.

# Tips, FAQ & Troubleshooting
- always first perform a dry-run with option -n
- if unsure why a certain rule will be executed use option --reason in the dry run, this will give the reason for the execution of each rule
- when there are 3 or less samples all the UMAP data and plots will be empty
- when there are 2 or less samples all the PCA data and plots will be empty
- in case the pipeline crashes, you manually canceled ll your jobs or when snakemake tries to "resume.. resubmit.." jobs, then remove .snakemake/incomplete directory!
- if you commit a lot of jobs eg via slurm (>500) this might take some time (ie 1s/job commit)
- two comments on peak support quantification
    - even though the peak support of a region in a certain sample is 0, does not mean that there are no reads counted in the count matrix, it just states that there was no peak called
    - the peak support can be >1 for certain samples in case of a consensus region that spans more than one region (with peaks) in the respective sample
- command for generating the workflow's rulegraph
```
snakemake --rulegraph --forceall | dot -Tsvg > workflow/dags/atacseq_pipeline_rulegraph.svg
```
provided in workflow/dags
- command for generating the directed acyclic graph (DAG) of all jobs with current configuration
```
snakemake --dag --forceall | dot -Tsvg > workflow/dags/all_DAG.svg
```
provided for both test examples in workflow/dags

# Appendix
## pipeline configuration specifications
- ...

## project configuration specifications
- ...
- 2 examples (hg38 & mm10) are given in the test/ directory

## sample annotation specifications
- every sequencing unit is a row and one sample can have multiple rows
- first column (sample_name) contains the sample name
- mandatory columns: data_source, skip_prepocess, all colums used in the respective data_source field in the pipeline configuration
- columns describe technical variables flowcell, lane,...
- 2 examples (hg38 & mm10) are given in the test/ directory


## metadata annotation specifications
- every sample is a row and the columns correspond to metadata entries eg cell type or condition
- mandatory columns: pass_qc with numeric value between 0 and 1; every sample >0 is included in the analysis
- 2 examples (hg38 & mm10) are given in the test/ directory


