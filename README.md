# Ultimate ATAC-seq Data Processing & Analysis Pipeline
From r**A**w (unaligned) BAM files to normali**Z**ed counts

A Snakemake implementation of the [BSF's](https://www.biomedical-sequencing.org/) [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline") extended by downstream processing and unsupervised analyses steps using bash, python and R. Reproducibility is ensured by using conda and Singularity.

![Workflow Rulegraph](./workflow/dags/atacseq_pipeline_rulegraph.svg)

# Features
- Processing
    - alignment of both single-end and paired-end reads in raw/unaligned BAM format (bowtie2)
    - peak calling (macs2)
    - peak annotation and motif analysis (homer)
    - MultiQC report generation (multiqc)
- Quantification
    - consensus region set generation
    - consensus region set annotation (uropa using regulatory build and gencode as refernce, and homer)
    - read count and peak support quantification of the consensus region set across samples, yielding a count and a support matrix with dimensions regions X samples
- optional: split data in multiple data subsets (eg by cell type, condition)
- Downstream Analysis (performed on the whole data set and each split separately)
    - region filtering
    - normalization by two different methods separately ([TMM](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25) & [CQN](https://academic.oup.com/biostatistics/article/13/2/204/1746212))
    - highly variable region (HVR) selection for each normalized data set
- Visualization after each analysis step
    - step specific plots (eg region filtering, HVR selection)
    - dimensionality reduction by PCA & UMAP and plotting of provided metadata
    - mean-variance-relationship plot
- Snakemake report generation for workflow statistics, documentation, reproducibility and result presentation

# Recommended Usage
1. perform only the processing, by setting the downstram_analysis flag to 0 in the project configuration
2. use the generated multiQC report (project_path/atacseq_report/multiqc_report.html) to judge the quality of your samples
3. fill out the mandatory quality control column (pass_qc) in the sample metadata configuration file (you can even use some of the quality metrics for plotting eg like I did in the example files with 'FRiP')
4. finally execute the remaining donwstream analysis steps by setting the respective flag to 1 only on the samples that passed quality control

# Installation (<10 minutes)
1. install snakemake, which requires conda & mamba, according to the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. clone/download this repository (eg git clone https://github.com/sreichl/atacseq_pipeline.git)

All software/package dependencies are installed and managed automatically via Snakemake and conda.

# Configuration
You need 2 configuration files and 2 annotation files to run the complete workflow from A to Z. You can use the provided examples as starting point. Always use absolute paths. If in doubt try the default values.
- pipeline configuration (conifg/pipeline_config.yaml): one off set up in your environment
    - before execution edit the path to your project configuration (project_config)
    - if you are a [CeMMie](https://cemm.at/) you can just take the provided file and only edit the project_config path
- project configuration (project_config): different for every project/dataset
- sample annotation (sample_annotation): technical unit annotation for processing
- sample metadata (sample_metadata): metadata used in downstream analysis steps and for visualization

detailed specifications can be found in the appendix below.

# Execution
## 1. Change working directory & activate conda environment
Execute always from within top level of the pipeline directory (ie atacseq_pipeline/).
Snakemake commands only work from within the snakemake conda environment.
```
cd atacseq_pipeline
conda activate snakemake
```
## 2. Execute a dry-run
command for a dry-run with option -n (-p makes Snakemake print the resulting shell command for illustration)
```
snakemake -p -n
```
## 3. Execute workflow local or on a cluster
### 3a. Local execution
command for execution with one core
```
snakemake -p -j1 --use-conda
```
### 3b. Cluster execution
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

## X. Singularity execution (not tested)
Singularity has to be installed (system wide by root) and available/loaded (eg module load singularity).
The pipeline automatically loads the correct singularity image from [Dockerhub](https://hub.docker.com/r/sreichl/atacseq_pipeline)

command for execution with singularity, just add the flag --use-singularity and use --singularity-args to provide all the necessary directories the pipeline needs access to (in the example it is configured for the three relevant partitions at CeMM)
```
snakemake -p --use-conda --use-singularity --singularity-args "--bind /nobackup:/nobackup --bind /research:/research --bind /home:/home"
```
# Use as module in another Snakemake workflow (soon)
- [https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules)
- [https://slides.com/johanneskoester/snakemake-6#/8](https://slides.com/johanneskoester/snakemake-6#/8)

# Report
command for report generation (this can take a few minutes, depending on the size of the dataset)
```
snakemake --report /absolute/path/to/report.zip
```

The command creates a self contained HTML based report in a zip archive containing the following sections:
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
- split1 (optional): Downstream analysis results of subset 1 of the whole data
- split2 (optional): Downstream analysis results of subset 2 of the whole data

# Examples
We provide configuration files for two example datasets (mm10 & hg38).
Additionally, the report zip archive of the hg38 test example is provided to showcase the pipeline results.

# Tips, FAQ & Troubleshooting
- always first perform a dry-run with option -n
- if unsure why a certain rule will be executed use option --reason in the dry run, this will give the reason for the execution of each rule
- when there are 3 or less samples all the UMAP data and plots will be empty
- when there are 2 or less samples all the PCA data and plots will be empty
- in case the pipeline crashes, you manually canceled your jobs or when snakemake tries to "resume.. resubmit.." jobs, then remove the .snakemake/incomplete directory!
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
The pipeline configuration has to be set up only once in your environment and can be found in config/pipeline_config.yaml. It is mostly home to a lot of absolute paths pointing to required reference data.
- project_config: absolute path to your project configuration file (see next section)
- adapter_fasta: (nextera) adapter fasta (.fa) file

the following files have always to be provided for each genome (hg39 & mm10), respectively

- bowtie2_index: indices for Bowtie2
- chromosome_sizes: files to define the chromosome lengths for a given genome
- blacklisted_regions: blacklisted regions from ENCODE as .bed files
- whitelisted_regions: complement to the blacklisted regions, also as .bed files
- unique_tss: .bed files, from gencode (hg38) or CCDS (mm10)
- regulatory_regions: regulatory build regulatory features chromosomes only from eg Ensembl
- genome_sizes: length of genomes as integers
- mitochondria_names: abbreviation of mitochondria chromosomes
- genome_fasta: genomes as .fa files
- gencode_gtf: gencode .gtf files
- regulatory_build_gtf: regulatory build regulatory features .gtf files (generated from gff.gz files from Ensembl with provided script /workflow/scripts/parse_reg_build_file.py)
- data_sources: here you can specify the dynamic absolute paths to your raw data by using {variables}, which have to be columns in the sample annotation configuration file. The entries correspond to the values in the data_source column in sample annotation configuration file.
- cluster parameters:
    - partition: on which partition to submit the jobs
    - memory: how much memory for each job
    - threads: how many threads/cpus-per-task are needed for each job

the provided configuration file is completely filled and set up for the environment at CeMM 

## project configuration specifications
- general project information
    - project_name: project name, used throughout the pipeline
    - project_uuid: unique ID
    - public_html_folder: in the folder a symlink with the unique ID will be created pointing to the results directory
    - base_url: url from where the public_html directory can be accessed
    - sample_annotation: path to the sample annotation file (see below)
    - project_path: path to the project directory for all the generated results
    - project_config: path to this file
    - genome: used genome in the project (hg38 or mm10)
    - adapter_sequence: nucleotide adapter sequence of the used ATAC-seq protocol (used by bowtie2 if provided, can be empty)
- downstream analysis parameters
    - downstream_analysis: 0 for stopping after processing (and inspection of multiQC report in project_path/atacseq_report/multiqc_report.html), 1 for executing the whole pipeline including all downstream steps
    - sample_metadata: path to sample metadata file (see below)
    - plot_by: list of metadata of interest to be plotted, elements need to be columns in the sample metadata file (see below) (put only 'pass_qc' for simple plots)
    - peak_support_threshold: minimum number of peaks per region, for region to be kept after filtering 
    - proportion: proportion of samples within min_group that has to "express" a region to not be filtered out
    - min_group: column in the sample metadata file that distinguishes between groups of interest (eg condition), used during region filtering step to ensure that regions of interest, which are only present in subgroups, are retained (can be empty, eg '')
    - split_by: column in the sample metadata file that distinguishes between subsets of the data that should be analyzed individually
    - HVR_percentage: percentage of regions to keep as highly variable regions
- region annotation parameters (number of bases, all integers)
    - tss_size: assumed size of transcription start sites (TSS)
    - proximal_size_up: assumed TSS proximal distance upstream
    - proximal_size_dn: assumed TSS proximal distance downstream
    - distal_size: distal distance 
- MultiQC report parameters (self explanatory, when in doubt go with default)
- UCSC Genome Browser trackhubs parameters (self explanatory, when in doubt go with default)

2 examples (hg38 & mm10) are provided in the test/ directory

## sample annotation specifications
- every sequencing unit is a row and one sample can have multiple rows
- first column (sample_name) contains the sample name
- mandatory columns: data_source, skip_prepocess (yes/no, if yes -> sample will not be processed), all columns used in the respective data_source field in the pipeline configuration
- columns describe technical variables flowcell, lane,...

2 examples (hg38 & mm10) are provided in the test/ directory


## sample metadata specifications
- every sample is a row and the columns correspond to metadata entries eg cell type or condition
- first column (sample_name) contains the sample name
- mandatory columns: pass_qc with numeric value between 0 and 1; every sample >0 is included in the downstream analysis

2 examples (hg38 & mm10) are provided in the test/ directory
