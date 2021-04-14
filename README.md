# atacseq_pipeline
Snakemake implementation of the BSF's [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline")


# snakemake commands

cmd for rulegraph
```
snakemake --rulegraph --forceall | dot -Tsvg > workflow/dags/atacseq_pipeline_rulegraph.svg
```

cmd for DAG with all jobs with current config
```
snakemake --dag --forceall | dot -Tsvg > workflow/dags/all_DAG.svg
```

cmd for cluster execution -> submits everything as separate job with dependencies
```
snakemake -p --profile config/slurm.cemm --use-conda --reason
```

flag for configfile 
```
--configfile path/to/project_configfile.yaml
```