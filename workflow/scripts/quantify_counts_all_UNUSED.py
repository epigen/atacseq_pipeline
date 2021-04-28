#!/bin/env python

import os
import pandas as pd
import pybedtools as bedtools

# configuration
output = snakemake.output[0]
chrom_file = snakemake.params["chrom_file"]
results_dir = snakemake.params["results_dir"]
# load annotation file and get sample names
annotations = pd.read_csv(snakemake.input[0], index_col=0)
annotations=annotations[(annotations['pass_qc']>0)]

elements_to_quantify = bedtools.BedTool(snakemake.input[1])

results=[]

for sample in annotations.index:
    bamfile=os.path.join(results_dir,"{}".format(sample),"mapped", "{}.filtered.bam".format(sample))
    print("Processing "+sample)
    try:
        result= elements_to_quantify.coverage(b=bamfile,sorted=True,g=chrom_file).to_dataframe(
                    names=["CHR", "START", "END", "ID", sample, "NA1", "NA2", "NA3"],
                    dtype={sample: int},
                    usecols=['ID', sample],
                    index_col='ID').T
        results.append(result)
    except Exception as e:
        print("Error occured while processing sample "+sample)
        results.append(pd.DataFrame(0,index=elements_to_quantify.index,columns=[sample]).T)

results = [item for item in results if item is not None]
results = pd.concat(results)
results.T.to_csv(output,index_label='ID')