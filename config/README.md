# Configuration

You need one configuration file and one annotation file to run the complete workflow. You can use the provided examples as starting point. If in doubt read the comments in the config and/or try the default values.

- project configuration (`config/config.yaml`): different for every project/dataset and configures the processing and quantification. The fields are described within the file.
- annotation.csv: CSV file consisting of one technical sequencing unit per row (i.e., one sample can include multiple sequencing units, hence mutliple rows) and 4 mandatory columns:
  - sample_name (first column!)
  - read_type: "single" or "paired".
  - bam_file: path to the raw/unaligned/unmapped [uBAM](https://gatk.broadinstitute.org/hc/en-us/articles/360035532132-uBAM-Unmapped-BAM-Format) file.
  - pass_qc: number between 0 (not used for downstream steps e.g., quantification) and 1. Every sample with pass_qc>0 is included in the downstream quantification and annotation steps.
  - (optional) additional sample metadata columns can be added and indicated for inclusion in the report.

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.

Two examples (hg38 & mm10) are provided in the `test/` directory and are pre-filled with resource locations from the [respective Zenodo download](../README.md#resources).
