#!/bin/env python

#### libraries
import os
import pandas as pd
import pybedtools as bedtools

#### configurations

# input
consensus_regions_path = snakemake.input["consensus_regions"]
bamfile_path = snakemake.input["bamfile"]
chrom_file = snakemake.config["chromosome_sizes"]

# output
quant_count_path = snakemake.output["quant_counts"]

# parameters
sample = snakemake.wildcards["sample"]

elements_to_quantify = bedtools.BedTool(consensus_regions_path)

print("Processing "+sample)
try:
    result = elements_to_quantify.coverage(b=bamfile_path,sorted=True,g=chrom_file).to_dataframe(
                names=["CHR", "START", "END", "ID", sample, "NA1", "NA2", "NA3"],
                dtype={sample: int},
                usecols=['ID', sample],
                index_col='ID').T
    result.to_csv(quant_count_path)
except Exception as e:
    print("Error occured while processing sample "+sample)
    elements_to_quantify_df = elements_to_quantify.to_dataframe(
                names=["CHR", "START", "END", "ID"],
                index_col='ID')
    pd.DataFrame(0,index=elements_to_quantify_df.index,columns=[sample]).T.to_csv(quant_count_path)