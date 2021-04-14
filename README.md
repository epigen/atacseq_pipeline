# atacseq_pipeline
Snakemake implementation of the BSF's [ATAC-seq Data Processing Pipeline](https://github.com/berguner/atacseq_pipeline "ATAC-seq Data Processing Pipeline")


# snakemake commands (always from within the snakemake conda environment!)

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
snakemake -p --profile config/slurm.cemm --use-conda --reason --configfile path/to/project_configfile.yaml
```

if more than 25 jobs (current limit of submitted jobs at the same time) then submit the command also as a job
```
sbatch -J atacseq_pipeline --mem=16000 --partition=longq --qos=longq --time=14-00:00:00 --wrap='snakemake -p --profile config/slurm.cemm --use-conda --reason --configfile path/to/project_configfile.yaml'
```
