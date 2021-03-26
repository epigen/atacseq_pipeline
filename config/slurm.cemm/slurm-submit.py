#!/usr/bin/env python3
"""
Snakemake SLURM submit script.
"""
import warnings  # use warnings.warn() rather than print() to output info in this script
import sys

from snakemake.utils import read_job_properties

import slurm_utils

# cookiecutter arguments
CLUSTER_CONFIG = "cluster.config.yaml"
ADVANCED_ARGUMENT_CONVERSION = {"yes": True, "no": False}[
    "no"
]

RESOURCE_MAPPING = {
    "time": ("time", "runtime", "walltime"),
    "mem": ("mem", "mem_mb", "ram", "memory"),
    "mem-per-cpu": ("mem-per-cpu", "mem_per_cpu", "mem_per_thread"),
    "nodes": ("nodes", "nnodes"),
}

# parse job
args = slurm_utils.parse_jobscript()

jobscript = args.jobscript
dependencies = args.dependencies

job_properties = read_job_properties(jobscript)

sbatch_options = {}
cluster_config = slurm_utils.load_cluster_config(CLUSTER_CONFIG)

# 2) cluster_config defaults
sbatch_options.update(cluster_config["__default__"])

# 3) Convert resources (no unit conversion!) and threads
sbatch_options.update(
    slurm_utils.convert_job_properties(job_properties, RESOURCE_MAPPING)
)

# 4) cluster_config for particular rule
sbatch_options.update(cluster_config.get(job_properties.get("rule"), {}))

# 5) cluster_config options
sbatch_options.update(job_properties.get("cluster", {}))

print("Dependencies are: {}".format(dependencies), file=sys.stderr)
dependency_option=""
if (dependencies is not None and len(dependencies)!=0):
    sbatch_options['dependency']="afterok:"+":".join(dependencies)

# 7) Format pattern in snakemake style
sbatch_options = slurm_utils.format_values(sbatch_options, job_properties)

# ensure sbatch output dirs exist
for o in ("output", "error"):
    slurm_utils.ensure_dirs_exist(sbatch_options[o]) if o in sbatch_options else None

# submit job and echo id back to Snakemake (must be the only stdout)
print(slurm_utils.submit_job(jobscript, **sbatch_options))
