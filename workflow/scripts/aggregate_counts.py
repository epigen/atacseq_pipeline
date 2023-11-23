#!/bin/env python

#### libraries
import os
import pandas as pd

#### configurations

# input
sample_paths = snakemake.input

# output
agg_path = snakemake.output[0]

# aggregate counts
results=[]

for sample_path in sample_paths:
    results.append(pd.read_csv(os.path.join(sample_path)))

results = pd.concat(results)
results.set_index(results.columns[0]).T.to_csv(agg_path)
