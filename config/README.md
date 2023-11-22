# Configuration

## NEW
required columns:
- sample_name (first column!)
- read_type: single or paired
- bam_file: path to bam file
- pass_qc: between 0 (not used for quantification) until 1 -> maybe use this as implicit quantification flag?

## OLD

You need 2 configuration files and 2 annotation files to run the complete workflow from A to Z. You can use the provided examples as starting point. Always use absolute paths. If in doubt try the default values.
- pipeline configuration (conifg/pipeline_config.yaml): one off set up in your environment
    - before execution edit the path to your project configuration (project_config)
    - if you are a [CeMMie](https://cemm.at/) you can just take the provided file and only edit the project_config path
- project configuration (project_config): different for every project/dataset
- sample annotation (sample_annotation): technical unit annotation for processing
- sample metadata (sample_metadata): metadata used in downstream analysis steps and for visualization


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
    - public_html_folder: in the folder a symlink with the unique ID will be created pointing to the results directory (this folder has to exist!)
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
- mandatory columns:
    - read_type (single/paired)
    - organism (mouse/human)
    - data_source: this field connects to the "data_sources" field in the pipeline configuration for dynamic absolute paths
    - skip_preprocess (yes/no, if yes -> sample will not be processed, this value has to be the same for all rows of the same sample!)
    - all columns used in the respective data_source field in the pipeline configuration
- columns describe technical variables flowcell, lane,...

2 examples (hg38 & mm10) are provided in the test/ directory


## sample metadata specifications
- every sample is a row and the columns correspond to metadata entries eg cell type or condition
- first column (sample_name) contains the sample name
- mandatory columns: pass_qc with numeric value between 0 and 1; every sample >0 is included in the downstream analysis


2 examples (hg38 & mm10) are provided in the test/ directory
