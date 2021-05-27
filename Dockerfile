FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="bb9a9b136a7a3c728c253049483eb5f9e7596dd0331682830a1416bc366f5a88"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/atacseq_analysis.yaml
#   prefix: /conda-envs/7446df8c60aaf3ad74baf0bb53bf99ff
#   name: atacseq_analysis
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.8
#     - bedtools==2.27.1
#     - ucsc-twobittofa
#     - bioconductor-deseq2
#     - bioconductor-cqn
#     - ipykernel
#     - bioconductor-copywriter
#     - rpy2
#     - sra-tools
#     - umap-learn
#     - pandas
#     - pybedtools
#     - seaborn
#     - matplotlib
#     - scikit-learn
#     - uropa
RUN mkdir -p /conda-envs/7446df8c60aaf3ad74baf0bb53bf99ff
COPY workflow/envs/atacseq_analysis.yaml /conda-envs/7446df8c60aaf3ad74baf0bb53bf99ff/environment.yaml

# Conda environment:
#   source: workflow/envs/atacseq_pipeline.yaml
#   prefix: /conda-envs/f2ab82afce323fed86ecb58f1eda1cb2
#   name: atacseq_pipeline
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - jmcmurray
#     - defaults
#   dependencies:
#     - samtools
#     - fastp
#     - samblaster=0.1.24
#     - bedtools
#     - macs2
#     - deeptools
#     - bowtie2
#     - setuptools
#     - json
#     - csvkit
#     - pandas
#     - matplotlib=3.2.2
RUN mkdir -p /conda-envs/f2ab82afce323fed86ecb58f1eda1cb2
COPY workflow/envs/atacseq_pipeline.yaml /conda-envs/f2ab82afce323fed86ecb58f1eda1cb2/environment.yaml

# Conda environment:
#   source: workflow/envs/multiqc.yaml
#   prefix: /conda-envs/e752bf04afb9ee98d1f1d66070f7e393
#   name: multiqc
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - multiqc=1.9
#     - setuptools
#     - pandas
#     - pip
#     - pip:
#       - ../../workflow/scripts/multiqc_atacseq
RUN mkdir -p /conda-envs/e752bf04afb9ee98d1f1d66070f7e393
COPY workflow/envs/multiqc.yaml /conda-envs/e752bf04afb9ee98d1f1d66070f7e393/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/7446df8c60aaf3ad74baf0bb53bf99ff --file /conda-envs/7446df8c60aaf3ad74baf0bb53bf99ff/environment.yaml && \
    mamba env create --prefix /conda-envs/f2ab82afce323fed86ecb58f1eda1cb2 --file /conda-envs/f2ab82afce323fed86ecb58f1eda1cb2/environment.yaml && \
    mamba env create --prefix /conda-envs/e752bf04afb9ee98d1f1d66070f7e393 --file /conda-envs/e752bf04afb9ee98d1f1d66070f7e393/environment.yaml && \
    mamba clean --all -y
