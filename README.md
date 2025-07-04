[![MrBiomics](https://img.shields.io/badge/MrBiomics-red)](https://github.com/epigen/MrBiomics/)
[![DOI](https://zenodo.org/badge/350342694.svg)](https://zenodo.org/doi/10.5281/zenodo.6323634)
[![](https://tokei.rs/b1/github/epigen/atacseq_pipeline?category=code)]() 
[![](https://tokei.rs/b1/github/epigen/atacseq_pipeline?category=files)]()
[![GitHub license](https://img.shields.io/github/license/epigen/atacseq_pipeline)](https://github.com/epigen/atacseq_pipeline/blob/master/LICENSE)
![GitHub Release](https://img.shields.io/github/v/release/epigen/atacseq_pipeline)
[![Snakemake](https://img.shields.io/badge/Snakemake->=8.20.1-green)](https://snakemake.readthedocs.io/en/stable/)

# Ultimate ATAC-seq Data Processing, Quantification & Annotation Pipeline
A [Snakemake 8](https://snakemake.readthedocs.io/en/stable/) workflow implementation of the [BSF's](https://www.biomedical-sequencing.org/) [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline") extended by downstream quantification and annotation steps using Bash and Python.

> [!NOTE]  
> This workflow adheres to the module specifications of [MrBiomics](https://github.com/epigen/MrBiomics), an effort to augment research by modularizing (biomedical) data science. For more details, instructions, and modules check out the project's repository.
>
> ⭐️ **Star and share modules you find valuable** 📤 - help others discover them, and guide our future work!

> [!IMPORTANT]  
> **If you use this workflow in a publication, please don't forget to give credit to the authors by citing it using this DOI [10.5281/zenodo.6323634](https://doi.org/10.5281/zenodo.6323634).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

# 🖋️ Authors
- [Stephan Reichl](https://github.com/sreichl)
- [Bekir Ergüner](https://github.com/berguner)
- [Daniele Barreca](https://github.com/DanieleBarreca)
- [Lukas Folkman](https://github.com/lukas-folkman)
- [Fangwen Zhao](https://github.com/fwzhao)
- [Rob ter Horst](https://github.com/rubbert)
- [Lina Dobnikar](https://github.com/ld401)
- [Christoph Bock](https://github.com/chrbock)

# 💿 Software
This project wouldn't be possible without the following software and their dependencies:

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

# 🔬 Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (`workflow/envs/*.yaml file`) or post-execution in the result directory (`atacseq_pipeline/envs/*.yaml`). Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g., [X].

**Processing.**
Sequencing adapters were removed using the software fastp (ver) [ref]. Bowtie2 (ver) [ref] was used for the alignment of the short reads (representing locations of transposition events) to the [GRCh38 (hg38)/GRCm38 (mm10)] assembly of the [human/mouse] genome using the “--very-sensitive” parameter. PCR duplicates were marked using samblaster (ver) [ref]. Aligned BAM files were then sorted, filtered using ENCODE blacklisted regions [ref], samtools view flags [SAM_flag], and indexed using samtools (ver) [ref]. To detect the open chromatin regions, peak calling was performed using MACS2 (ver) [ref] using the “--nomodel”, “--keep-dup [macs2_keep_dup]” and “--extsize 147” options on each sample. HOMER (ver) [ref] function findMotifs was used for motif enrichment analysis of the detected open chromatin regions. Quality control metrics were aggregated and reported using MultiQC (ver) [ref], [X] sample(s) needed to be removed.

**Quantification.**
A consensus region set, comprising of [X] genomic regions, was generated, by merging the identified peak summits, extended by [slop_extension]bp on both sides using the slop function from bedtools (ver) [ref] and pybedtools (ver) [ref], across all samples while again discarding peaks overlapping blacklisted features as defined by the ENCODE project [ref]. The consensus region set was used to quantify the chromatin accessibility in each sample by summing the number of reads overlapping each consensus region. The consensus region set, and sample-wise quantification of accessibility was performed using bedtools (ver) [ref] and pybedtools (ver) [ref]. Furthermore, the consensus region set was used to quantify the peak support per sample and each region was mapped to its closest TSS according to the HOMER annotation within proximal TSS up and downstream distances [proximal_size_up/down] yielding a gene/TSS-based quantification. Complementary, all promoter regions, defined by the same parameters, were quantified for each sample and aggregated to yield a gene/promoter-based quantification. Finally, all sample-wise enriched known motifs according to HOMER were aggregated.

**Annotation.**
Consensus regions were annotated using annotatePeaks function from HOMER (ver) [ref]. Additionally, we annotated all consensus regions using UROPA (ver) [ref] and genomic features from the [GENCODE vX] basic gene annotation as: “TSS proximal” if the region’s midpoint was within [X] bp upstream or [X] bp downstream from a TSS, or if the region overlapped with a TSS; “gene body” if the region overlapped with a gene; “distal” if the region’s midpoint was within [X] bp of a TSS; and “intergenic” otherwise. For each region, only the closest feature was considered, and the annotations took precedence in the following order: TSS proximal, gene body, distal, and intergenic. Finally, bedtools was employed to quantify nucleotide counts and proportional content per consensus region.

The processing and quantification described here was performed using a publicly available Snakemake [ver] (ref) workflow [[10.5281/zenodo.6323634](https://doi.org/10.5281/zenodo.6323634)].

# 🚀 Features
- Processing (`results/`)
    - Alignment of both single-end and paired-end reads in raw/unaligned/unmapped [uBAM](https://gatk.broadinstitute.org/hc/en-us/articles/360035532132-uBAM-Unmapped-BAM-Format) format with Bowtie2.
      - Filtering using `samtools view` can be configured using [SAM Flags](https://broadinstitute.github.io/picard/explain-flags.html) (`SAM_flag`).
    - Peak calling with `MACS2`.
      - Duplicate handling can be configured using `macs2_keep_dup`.
      - Even though the peak support of a region in a certain sample is 0, does not mean that there are no reads counted in the count matrix, it just states that there was no peak called.
      - The peak support can be >1 for certain samples in case of a consensus region spanning more than one peak within the respective sample.
    - Peak annotation and motif analysis HOMER.
    - Quantification of TSS coverage.
- Reporting (`report/`)
    - MultiQC report generation using MultiQC, extended with an in-house developed plugin [atacseq_report](./workflow/scripts/multiqc_atacseq).
- Quantification (`counts/`)
    - Consensus region set generation across all called peaks (`consensus_regions.bed`).
    - Read count quantification of the consensus regions across samples, yielding a count matrix with dimensions consensus regions X samples (`consensus_counts.csv`).
    - Peak support quantification of the consensus regions across samples, yielding a count matrix with dimensions consensus regions X samples (`support_counts.csv`).
    - Consensus regions mapped to closest gene TSS according to HOMER (Distance to TSS) within proximal TSS up and downstream distances (`TSS_regions.bed`, `TSS_counts.csv`, `TSS_annotation.csv`).
    - Read count quantification of promoter regions based on provided proximal TSS up and downstream distances (`promoter_regions.bed`, `promoter_counts.csv`, `promoter_annotation.csv`).
      - [Pseudoautosomal regions in human](https://www.ensembl.org/info/genome/genebuild/human_PARS.html) chromosome `Y` are skipped.
    - Aggregation of all sample-wise HOMER known motif enrichment results into one CSV in long-format (`HOMER_knownMotifs.csv`).
- Annotation (`counts/`)
    - Sample annotation file based on MultiQC general stats and provided annotations for downstream analysis (`sample_annotation.csv`).
    - Consensus region set annotation using (`consensus_annotation.csv`)
      - `UROPA` with regulatory build and gencode as references, configurable here: `workflow/resources/UROPA/*.txt`.
      - `HOMER` with `annotatePeaks.pl`. NB: We have empirically found, that some human sex genes, e.g., the well established protein coding genes UTY and STS, are not annotated.
      - `bedtools` for nucleotide counts/content (e.g., % of GC).

> [!IMPORTANT]  
> **Duplciate reads** can be filtered during the alignment step by `samtools` and/or ignored during peak calling by `MACS2`.
> **The inclusion of duplicates** should be intentional, and may lead to a large number of consensus regions.
> **The removal of duplicates** should be intentional, might remove real biological signal.
> **The decision depends** on your downstream analysis steps e.g., rigorous filtering (e.g., using `edgeR::filterByExpr`) and/or accounting for PCR bias by normalization conditional on genomic region length and GC content (e.g., [CQN](https://academic.oup.com/biostatistics/article/13/2/204/1746212)) and goals (e.g., differential accessibility analysis).
> **We recommend** reading this ChIP-seq tutorial's section on ["Removing redundancy"](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html).

# 🛠️ Usage
These steps are the recommended usage for this workflow:

0. Configure the workflow by pointing to the relevant resources, e.g., downloaded from Zenodo for [hg38 or mm10 (see instructions below)](#resources).
1. Perform only the processing, by setting the pass_qc annotation for all samples to 0.
2. Use the generated MultiQC report (result_path/ataceq_pipeline/report/multiqc_report.html) to judge the quality of each sample (see tips in the next section on [Quality Control](#quality-control)).
3. Fill out the mandatory quality control column (pass_qc) in the annotation file accordingly (everything >0 will be included in the downstream steps).
4. Finally, execute the remaining downstream quantification and annotation steps by running the workflow. Thereby only the samples that passed quality control will be included in the consensus region set generation (i.e., the feature space) and all downstream steps.

This workflow is written with Snakemake and its usage is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/atacseq_pipeline).

# ⚙️ Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# 📖 Examples
Explore a detailed example showcasing module usage and downstream analysis in our comprehensive end-to-end [MrBiomics Recipe](https://github.com/epigen/MrBiomics?tab=readme-ov-file#-recipes) for [ATACseq Analysis](https://github.com/epigen/MrBiomics/wiki/ATAC%E2%80%90seq-Analysis-Recipe), including data, configuration, annotation and results.

# 🔍 Quality Control
Below are some guidelines for the manual quality control of each sample, but keep in mind that every experiment/dataset is different.

1. Reads Mapped ~ $30\cdot 10^{6}$ ($>20\cdot 10^{6}$ at least)
2. % Aligned >90%
3. % Mitochondrial <10%
4. Peaks (depend on reads)
    - FriP (Fraction of reads in Peaks) ~ >20% (can be misleading as 80-90% are also not good)
    - Regulatory regions >10% (as it is roughly 10% of the genome)
    - TSS (Transcription Start Site) normalized coverage ideally > 4 (at least >2)
    - % Duplications “not excessive”
5. Inspect [Genome Browser Tracks](https://github.com/epigen/genome_tracks/) using UCSC Genome Browser (online) or IGV (local)
    - Compare all samples to the best, based on above's QC metrics.
    - Check cell type / experiment-specific markers or sex chromosome (`X`/`Y`) for accessibility as positive controls.
    - Check e.g., developmental regions for accessibility as negative controls.
6. [Unsupervised Analysis](https://github.com/epigen/unsupervised_analysis) (e.g., PCA or UMAP)
    - Identify outliers/drivers of variation, especially in the control samples and within replicates.

> [!IMPORTANT]  
> Sometimes reads map to `Y` in females, because `X` and `Y` chromosomes both have [pseudoautosomal regions (PARs)](https://www.ensembl.org/info/genome/genebuild/human_PARS.html) that are common between the two chromosomes.
  
My personal QC value scheme to inform downstream analyses (e.g., unsupervised analysis)
- 0 = did not pass
- 2 options
  - for every metric that is not ok subtract 0.25 from 1, which means it requires 4 “strikes” for a sample to be removed due to low quality.
  - alternative
      - 0.5 = passed with reservations (e.g., metrics and genome browser tracks were not optimal, but still good enough)
      - 0.75 = not ideal (e.g., at least metrics or IGV tracks were not optimal)
- 1 = passed (perfect)

Finally, a previous PhD student in our lab, [André Rendeiro](https://orcid.org/0000-0001-9362-5373), wrote about ["ATAC-seq sample quality, wet lab troubleshooting and advice"](https://github.com/epigen/open_pipelines/blob/master/pipelines/atacseq.md#sample-quality-wet-lab-troubleshooting-and-advice).

# 🔗 Links
- [GitHub Repository](https://github.com/epigen/atacseq_pipeline/)
- [GitHub Page](https://epigen.github.io/atacseq_pipeline/)
- [Zenodo Repository](https://doi.org/10.5281/zenodo.6323634)
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/atacseq_pipeline)

# 📚 Resources
- Data Resources: To ensure the reproducibility of results and to make the workflow accessible we provide all required reference data for the analysis of ATAC-seq samples for [human GRCh38 (hg38)](https://doi.org/10.5281/zenodo.6344173) and [mouse GRCm38 (mm10)](https://doi.org/10.5281/zenodo.6344321) genomes on Zendodo.

  **Command line**
  ```console
  # download Zenodo records using zenodo_get

  # install zenodo_get v1.3.4
  conda install -c conda-forge zenodo_get=1.3.4

  # human GRCh38 (hg38)
  zenodo_get --record 6344173 --output-dir=resources/atacseq_pipeline/hg38/
  cd resources/atacseq_pipeline/hg38
  unzip indices_for_Bowtie2.zip && rm indices_for_Bowtie2.zip

  # mouse GRCm38 (mm10)
  zenodo_get --record 6344322 --output-dir=resources/atacseq_pipeline/mm10/
  cd resources/atacseq_pipeline/mm10
  unzip indices_for_Bowtie2.zip && rm indices_for_Bowtie2.zip
  ```
  **Snakemake rule** for workflows
  ```python
  #### Get resources from Zenodo as custom Snakemake rule ####
  # Downloads Bowtie2 indices for hg38 from Zenodo record 6344173 and unpacks them.
    rule MyATAC_get_resources:
        output:
            "resources/MyATAC/atacseq_pipeline/hg38/gencode.v38.basic.annotation.gtf",
            resource_dir = directory("resources/MyATAC/atacseq_pipeline/hg38/"),
        params:
            zenodo_record = "6344173",
            zip_filename = "indices_for_Bowtie2.zip"
        conda:
            "../envs/zenodo_get.yaml"
        shell:
            """
            # Download the specific record to the target directory
            zenodo_get --record {params.zenodo_record} --output-dir={output.resource_dir}
    
            # Change directory, unzip the specific file, and remove the zip archive
            # Using && ensures commands run sequentially and stop if one fails
            cd {output.resource_dir} && \
            unzip {params.zip_filename} && \
            rm {params.zip_filename}
            """
  ```
- Recommended compatible [MrBiomics](https://github.com/epigen/MrBiomics) modules for
  - upstream analysis:
      - [Fetch Public Sequencing Data and Metadata Using iSeq](https://github.com/epigen/fetch_ngs/) to retrieve and prepare public ATAC-ses data for downstream processing.
  - downstream analysis (in that order):
      - [Genome Browser Track Visualization](https://github.com/epigen/genome_tracks/) for quality control and visual inspection/analysis of genomic regions/genes of interest or top hits.
      - [<ins>Sp</ins>lit, F<ins>ilter</ins>, Norma<ins>lize</ins> and <ins>Integrate</ins> Sequencing Data](https://github.com/epigen/spilterlize_integrate/) after count quantification.
      - [Unsupervised Analysis](https://github.com/epigen/unsupervised_analysis) to understand and visualize similarities and variations between cells/samples, including dimensionality reduction and cluster analysis. Useful for all tabular data including single-cell and bulk sequencing data.
      - [Differential Analysis with limma](https://github.com/epigen/dea_limma) to identify and visualize statistically significantly different features (e.g., genes or genomic regions) between sample groups.
      - [Enrichment Analysis](https://github.com/epigen/enrichment_analysis) for biomedical interpretation of (differential) analysis results using prior knowledge.
    - [Introduction to ChIP-seq using high performance computing](https://hbctraining.github.io/Intro-to-ChIPseq/)

# 📑 Publications
The following publications successfully used this module for their analyses.
- [Casteels et al. (2022) Cell Reports - SMNDC1 links chromatin remodeling and splicing to regulate pancreatic hormone expression.](https://doi.org/10.1016/j.celrep.2022.111288)
- ...

# ⭐ Star History

[![Star History Chart](https://api.star-history.com/svg?repos=epigen/atacseq_pipeline&type=Date)](https://star-history.com/#epigen/atacseq_pipeline&Date)
