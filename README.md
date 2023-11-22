# Ultimate ATAC-seq Data Processing & Quantification Pipeline
[![DOI](https://zenodo.org/badge/350342694.svg)](https://zenodo.org/badge/latestdoi/350342694)

From r**A**w (unaligned) BAM files to count**Z**.
A Snakemake implementation of the [BSF's](https://www.biomedical-sequencing.org/) [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline") extended by downstream quantification and annotation steps using bash, python and R. Reproducibility and scalability is ensured by using Snakemake and conda.

This workflow adheres to the module specifications of [MR.PARETO](https://github.com/epigen/mr.pareto), an effort to augment research by modularizing (biomedical) data science. For more details and modules check out the project's repository.

**If you use this workflow in a publication, please don't forget to give credits to the authors by citing it using this DOI [10.5281/zenodo.6323634](https://doi.org/10.5281/zenodo.6323634).**

![Workflow Rulegraph](./workflow/dags/atacseq_pipeline_rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [Links](#links)
  * [Resources](#resources)
  * [Publications](#publications)

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
| deeptools      | https://doi.org/10.1093/nar/gkw257                |
| ENCODE         | https://doi.org/10.1038/s41598-019-45839-z        |
| fastp          | https://doi.org/10.1093/bioinformatics/bty560     |
| HOMER          | https://doi.org/10.1016/j.molcel.2010.05.004      |
| MACS2          | https://doi.org/10.1186/gb-2008-9-9-r137          |
| MultiQC        | https://doi.org/10.1093/bioinformatics/btw354     |
| pybedtools     | https://doi.org/10.1093/bioinformatics/btr539     |
| pandas         | https://doi.org/10.5281/zenodo.3509134            |
| samblaster     | https://doi.org/10.1093/bioinformatics/btu314     |
| samtools       | https://doi.org/10.1093/bioinformatics/btp352     |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |
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

The processing and analysis described here was performed using a publicly available Snakemake [ver] (ref) workflow [[10.5281/zenodo.6323634](https://doi.org/10.5281/zenodo.6323634)].

# Features
- Processing
    - alignment of both single-end and paired-end reads in raw/unaligned BAM format (bowtie2)
    - peak calling (macs2)
      - Note: even though the peak support of a region in a certain sample is 0, does not mean that there are no reads counted in the count matrix, it just states that there was no peak called
      - Note: the peak support can be >1 for certain samples in case of a consensus region that spans more than one region (with peaks) in the respective sample
    - peak annotation and motif analysis (homer)
    - MultiQC report generation (multiqc)
- Quantification
    - consensus region set generation
    - consensus region set annotation (uropa using regulatory build and gencode as refernce, and homer)
    - read count and peak support quantification of the consensus region set across samples, yielding a count and a support matrix with dimensions regions X samples
- Snakemake report generation for workflow statistics, documentation, reproducibility and result presentation using the same structure as the result files.

# Results (TODO: integrate into Features)
Project directory structure:
- all: Downstream analysis results of the whole data, including consensus region set and region annotation
- atacseq_hub: genome browser track files (.bigWig) for each sample
- atacseq_report: statistics and metrics from the processing part as input for the MultiQC report
- atacseq_results: one directory per sample with all the processing and quantification results
- projectName_report.zip: self contained HTML Snakemake report
- split1 (optional): Downstream analysis results of subset 1 of the whole data
- split2 (optional): Downstream analysis results of subset 2 of the whole data

# Usage
These steps are the recommended usage for this workflow:
1. Perform only the processing, by setting the quantification flag to 0 in the configuration file.
2. Use the generated multiQC report (project_path/atacseq_report/multiqc_report.html) to judge the quality of your samples.
3. Fill out the mandatory quality control column (pass_qc) in the sample annotation file.
4. Finally, execute the remaining donwstream quantification and annotation steps by setting the quantification flag to 1, thereby only the samples that passed quality control will be included in the consensus region set generation (i.e., the feature space).

This workflow is written with Snakemake and its usage is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/atacseq_pipeline).

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
We provide configuration files for two example datasets (mm10 & hg38).
Additionally, the report zip archive of the hg38 test example is provided to showcase the pipeline results.

# Links
- [GitHub Repository](https://github.com/epigen/atacseq_pipeline/)
- [GitHub Page](https://epigen.github.io/atacseq_pipeline/)
- [Zenodo Repository](https://doi.org/10.5281/zenodo.6323634)
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/atacseq_pipeline)

# Resources
- Data Resources: To ensure reproducibility of results and to make the pipeline easy-to-use we provide all required reference data for the analysis of ATAC-seq samples for both supported genomes on Zendodo: 
  - [resources for the GRCm38 (mm10) assembly of the mouse genome](https://doi.org/10.5281/zenodo.6344321)
  - [resources for the GRCh38 (hg38) assembly of the human genome](https://doi.org/10.5281/zenodo.6344173)
- Recommended [MR.PARETO](https://github.com/epigen/mr.pareto) modules for downstream analyses (in that order):
  - [<ins>Sp</ins>lit, F<ins>ilter</ins>, Norma<ins>lize</ins> and <ins>Integrate</ins> Sequencing Data](https://github.com/epigen/spilterlize_integrate/)
  - [Unsupervised Analysis](https://github.com/epigen/unsupervised_analysis) to understand and visualize similarities and variations between samples.
  - [Differential Analysis with limma](https://github.com/epigen/dea_limma) to identify and visualize statistically significant features between sample groups.
  - [Enrichment Analysis](https://github.com/epigen/enrichment_analysis)  for biomedical interpretation of results.
  - [Genome Tracks](https://github.com/epigen/genome_tracks) for visualization of genmoics regions of interest or top hits.

# Publications
The following publications successfully used this module for their analyses.
- [Casteels et al. (2022) Cell Reports - SMNDC1 links chromatin remodeling and splicing to regulate pancreatic hormone expression.](https://doi.org/10.1016/j.celrep.2022.111288)
- ...
