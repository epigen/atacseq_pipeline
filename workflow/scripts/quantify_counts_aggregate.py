#!/bin/env python

import os
import pandas as pd

# configuration
output = snakemake.output[0]
results_dir = snakemake.params["results_dir"]

# load annotation file and get sample names
annotations = pd.read_csv(snakemake.input[0], index_col=0)

results=[]

for sample in annotations.index:
    print(sample)
    result=pd.read_csv(os.path.join(results_dir,"{}".format(sample),"mapped", "{}_quantification_counts.csv".format(sample)))
    results.append(result)

results = [item for item in results if item is not None]
results = pd.concat(results)
results.set_index(results.columns[0]).T.to_csv(output)