#!/bin/env python

#### libraries
import os
import pandas as pd

#### configurations

# input
support_sample_paths = snakemake.input["quant_support"]

# output
support_agg_path = snakemake.output["quant_support_agg"]

# output = snakemake.output[0]
# results_dir = snakemake.params["results_dir"]

# # load annotation file and get sample names
# annotations = pd.read_csv(snakemake.input[0], index_col=0)

results=[]

# for sample in annotations.index:
#     print(sample)
#     result=pd.read_csv(os.path.join(results_dir,"{}".format(sample),"peaks", "{}_quantification_support.csv".format(sample)))
#     results.append(result)
    
for sample_path in support_sample_paths:
    results.append(pd.read_csv(os.path.join(sample_path)))

# results = [item for item in results if item is not None]
results = pd.concat(results)
results.set_index(results.columns[0]).T.to_csv(support_agg_path)
