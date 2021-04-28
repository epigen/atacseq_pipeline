# atacseq_pipeline
Snakemake implementation of the BSF's [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline")

execution currently always from within root level ie atacseq_pipeline folder (to be adapted)

# snakemake commands (always from within the snakemake conda environment!)

cmd for rulegraph
```
snakemake --rulegraph --forceall | dot -Tsvg > workflow/dags/atacseq_pipeline_rulegraph.svg
```

cmd for DAG with all jobs with current config
```
snakemake --dag --forceall | dot -Tsvg > workflow/dags/all_DAG.svg
```

cmd for execution with one core
```
snakemake -p -j1 --use-conda --reason --config project_config=path/to/project_configfile.yaml
```

cmd for cluster execution -> submits everything as separate job with dependencies
```
snakemake -p --profile config/slurm.cemm --use-conda --reason --config project_config=path/to/project_configfile.yaml
```
