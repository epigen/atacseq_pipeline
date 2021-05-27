#!/bin/env python

import os
import pandas as pd
import pybedtools as bedtools

# configuration
output = snakemake.output[0]
chrom_file = snakemake.params["chrom_file"]
results_dir = snakemake.params["results_dir"]
sample=snakemake.wildcards["sample"]

elements_to_quantify = bedtools.BedTool(snakemake.input[0])

bamfile=snakemake.input[1]

print("Processing "+sample)
try:
    result= elements_to_quantify.coverage(b=bamfile,sorted=True,g=chrom_file).to_dataframe(
                names=["CHR", "START", "END", "ID", sample, "NA1", "NA2", "NA3"],
                dtype={sample: int},
                usecols=['ID', sample],
                index_col='ID').T
    result.to_csv(output)
except Exception as e:
    print("Error occured while processing sample "+sample)
    pd.DataFrame(0,index=elements_to_quantify.index,columns=[sample]).T.to_csv(output)