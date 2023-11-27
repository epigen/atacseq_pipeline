#!/bin/env python

#### libraries
import os
import pandas as pd

#### configurations

# input
support_sample_paths = snakemake.input["quant_support"]

# output
support_agg_path = snakemake.output["quant_support_agg"]

results=[]

for sample_path in support_sample_paths:
    results.append(pd.read_csv(os.path.join(sample_path)))

results = pd.concat(results)
results.set_index(results.columns[0]).T.to_csv(support_agg_path)
