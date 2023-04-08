# Ultimate ATAC-seq Data Processing & Analysis Pipeline
[![DOI](https://zenodo.org/badge/350342694.svg)](https://zenodo.org/badge/latestdoi/350342694)

From r**A**w (unaligned) BAM files to normali**Z**ed counts.

A Snakemake implementation of the [BSF's](https://www.biomedical-sequencing.org/) [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline") extended by downstream processing and unsupervised analyses steps using bash, python and R. Reproducibility and Scalability is ensured by using Snakemake, conda and Singularity.

**If you use this workflow in a publication, please don't forget to give credits to the authors by citing it using this DOI [10.5281/zenodo.6323634](https://doi.org/10.5281/zenodo.6323634).**

![Workflow Rulegraph](./workflow/dags/atacseq_pipeline_rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Installation](#installation)
  * [Configuration](#configuration)
  * [Execution](#execution)
  * [Module](#module)
  * [Report](#report)
  * [Results](#results)
  * [Examples](#examples)
  * [Data Resources](#resources)
  * [Tips & FAQs](#tips)
  * [Links](#links)

# Authors
- [Stephan Reichl](https://github.com/sreichl)
- [Bekir Ergüner](https://github.com/berguner)
- [Daniele Barreca](https://github.com/DanieleBarreca)
- [Lukas Folkman](https://github.com/lukas-folkman)
- [Lina Dobnikar](https://github.com/ld401)
- [Christoph Bock](https://github.com/chrbock)

# Software
This project wouldn't be possible without the following software

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| bedtools       | https://doi.org/10.1093/bioinformatics/btq033     |
| Bowtie2        | https://doi.org/10.1038/nmeth.1923                |
| CQN            | https://doi.org/10.1093/biostatistics/kxr054      |
| deeptools      | https://doi.org/10.1093/nar/gkw257                |
| ENCODE         | https://doi.org/10.1038/s41598-019-45839-z        |
| fastp          | https://doi.org/10.1093/bioinformatics/bty560     |
| HOMER          | https://doi.org/10.1016/j.molcel.2010.05.004      |
| MACS2          | https://doi.org/10.1186/gb-2008-9-9-r137          |
| matplotlib     | https://doi.org/10.1109/MCSE.2007.55              |
| MultiQC        | https://doi.org/10.1093/bioinformatics/btw354     |
| pybedtools     | https://doi.org/10.1093/bioinformatics/btr539     |
| pandas         | https://doi.org/10.5281/zenodo.3509134            |
| samblaster     | https://doi.org/10.1093/bioinformatics/btu314     |
| samtools       | https://doi.org/10.1093/bioinformatics/btp352     |
| SCANPY         | https://doi.org/10.1186/s13059-017-1382-0.        |
| scikit-learn   | http://jmlr.org/papers/v12/pedregosa11a.html      |
| seaborn        | https://doi.org/10.21105/joss.03021               |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |
| TMM            | https://doi.org/10.1186/gb-2010-11-3-r25          |
| UMAP           | https://doi.org/10.21105/joss.00861               |
| UROPA          | https://doi.org/10.1038/s41598-017-02464-y        |

# Methods
This is a template for the **Methods** section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

**Processing.**
Sequencing adapters were removed using the software fastp (ver) [ref]. Bowtie2 (ver) [ref] was used for the alignment of the short reads (representing locations of transposition events) to the [GRCh38 (hg38)/GRCm38 (mm10)] assembly of the [human/mouse] genome using the “--very-sensitive” parameter. PCR duplicates were marked using samblaster (ver) [ref]. Aligned BAM files were then sorted, filtered using ENCODE blacklisted regions [ref], and indexed using samtools (ver) [ref]. To detect the open chromatin regions, peak calling was performed using MACS2 (ver) [ref] using the “--nomodel”, “--keep-dup auto” and “--extsize 147” options on each sample. Homer (ver) [ref] function findMotifs was used for motif enrichment analysis over the detected open chromatin regions.

Quality control metrics were aggregated and reported using MultiQC (ver) [ref], [X] sample(s) needed to be removed.

**Quantification.**
A consensus region set, comprising of [X] genomic regions, was generated, by merging the identified peak summits, extended by 250 bp on both sides using the slop function from bedtools (ver) [ref] and pybedtools (ver) [ref], across all samples while again discarding peaks overlapping blacklisted features as defined by the ENCODE project [ref].
Consensus regions were annotated using annotatePeaks function from Homer (ver) [ref].

Additionally, we annotated all consensus regions using UROPA (ver) [ref] and genomic features from the [GENCODE vX] basic gene annotation as: “TSS proximal” if the region’s midpoint was within [X] bp upstream or [X] bp downstream from a TSS, or if the region overlapped with a TSS; “gene body” if the region overlapped with a gene; “distal” if the region’s midpoint was within [X] bp of a TSS; and “intergenic” otherwise. For each region, only the closest feature was considered, and the annotations took precedence in the following order: TSS proximal, gene body, distal, and intergenic. 

The consensus region set was used to quantify the chromatin accessibility in each sample by summing the number of reads overlapping each consensus region. The consensus region set, and sample-wise quantification of accessibility was performed using bedtools (ver) [ref] and pybedtools (ver) [ref].

**Optional.** We split up of data into subsets according to the annotation [X] and performed all downstream analyses for each subset separately.

**Downstream Analysis.**
For all downstream analyses, we filtered the [X] consensus regions to [X] regions which had reads in at least [X] samples, were covered by at least [X] reads (normalized by median library size) in at least [X] proportion of samples of the smallest subsample group with potential signal determined with the annotation [X], and by at least [X] total reads across all samples.

Next, we determined highly variable regions (HVR) using an adaption of a SCANPY (ver) [ref] function highly_variable_genes with flavor 'cellranger', but instead of dispersion=variation/mean we use dispersion=standard deviation, because in ATAC-seq data the number of regions might be very large. This could lead to log(cpm) values and log(cpm)-means being negative, resulting in negative dispersion values which are meaningless (additionally avoiding division by 0 problems). Therefore we only employ the binning strategy by mean for stabilization, but not a "coefficient of variation" strategy.

Conditional quantile normalization cqn (ver) [ref] was performed on the filtered accessibility matrix using the region length and GC-content of the consensus regions as conditions, quantified using bedtools (ver) [ref] and pybedtools (ver) [ref].

Trimmed mean of M-values (TMM) normalization (ver) [ref] was performed on the filtered accessibility matrix.

**Unsupervised Analysis & Visualization.**
We applied both linear and non-linear unsupervised dimensionality reduction methods to normalized data to visualize emerging sample-wise patterns in two dimensions. We used Principal Component Analysis (PCA) [ref] from scikit-learn (ver) [ref] as the linear approach and Uniform Manifold Approximation projection (UMAP) from umap-learn (ver) [ref] with the correlation metric as the non-linear approach. For visualization matplotlib (ver) [ref] was used.

The processing and analysis described here was performed using a publicly available Snakemake [ver] (ref) workflow [[10.5281/zenodo.6323634](https://doi.org/10.5281/zenodo.6323634)].

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

# Usage
These steps are the recommended usage for the first run of the workflow:
1. perform only the processing, by setting the downstram_analysis flag to 0 in the project configuration
2. use the generated multiQC report (project_path/atacseq_report/multiqc_report.html) to judge the quality of your samples
3. fill out the mandatory quality control column (pass_qc) in the sample metadata configuration file (you can even use some of the quality metrics for plotting eg like I did in the example files with 'FRiP')
4. finally execute the remaining donwstream analysis steps by setting the respective flag to 1, thereby only the samples that passed quality control will be included

This workflow is written with snakemake and its usage is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/atacseq_pipeline).

# Installation
This should take less than 10 minutes.
1. install snakemake, which requires conda & mamba, according to the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. clone/download this repository (eg git clone https://github.com/epigen/atacseq_pipeline.git)

All software/package dependencies are installed and managed automatically via Snakemake and conda (and optional Singularity) and installed upon the first run of the workflow.

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

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
The pipeline automatically loads the correct singularity image from [Dockerhub](https://hub.docker.com/r/
/atacseq_pipeline)

command for execution with singularity, just add the flag --use-singularity and use --singularity-args to provide all the necessary directories the pipeline needs access to (in the example it is configured for the three relevant partitions at CeMM)
```
snakemake -p --use-conda --use-singularity --singularity-args "--bind /nobackup:/nobackup --bind /research:/research --bind /home:/home"
```
# Module
Use as module in another Snakemake workflow (soon)
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

# Resources
To ensure reproducibility of results and to make the pipeline easy-to-use we provide all required reference data for the analysis of ATAC-seq samples for both supported genomes on Zendodo: 
- [resources for the GRCm38 (mm10) assembly of the mouse genome](https://doi.org/10.5281/zenodo.6344321)
- [resources for the GRCh38 (hg38) assembly of the human genome](https://doi.org/10.5281/zenodo.6344173)

# Tips
Here are some tips for troubleshooting & FAQs:
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

# Links
- [GitHub Repository](https://github.com/epigen/atacseq_pipeline/)
- [GitHub Page](https://epigen.github.io/atacseq_pipeline/)
- [Zenodo Repository](https://doi.org/10.5281/zenodo.6323634)
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/atacseq_pipeline)

# Publications
The following publications successfully used this module for their analyses.
- [Casteels et al (2022) Cell Reports - SMNDC1 links chromatin remodeling and splicing to regulate pancreatic hormone expression.](https://doi.org/10.1016/j.celrep.2022.111288)
