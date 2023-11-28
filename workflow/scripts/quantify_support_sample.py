#!/bin/env python

#### libraries
import os
import pandas as pd
import pybedtools as bedtools

#### configurations

# input
consensus_regions_path = snakemake.input["consensus_regions"]
peakfile_path = snakemake.input["peakfile"]
chrom_file = snakemake.config["chromosome_sizes"]

# output
quant_support_path = snakemake.output["quant_support"]

# parameters
sample = snakemake.wildcards["sample"]

# quantify peak support within sample
consensus_peaks = bedtools.BedTool(consensus_regions_path)
consensus_peaks_df = bedtools.BedTool(consensus_regions_path).to_dataframe().set_index('name')

# peakfile=snakemake.input[1]
result = pd.DataFrame(0,index=consensus_peaks_df.index,columns=[sample])
    
try:
    sample_peaks = bedtools.BedTool(peakfile_path)
    result = consensus_peaks.intersect(
            sample_peaks,
            g=chrom_file, 
            wa=True,
            c=True
        ).to_dataframe(
            index_col='name',
            usecols=[3,4],
            names=['name',sample]
        )
except Exception as e:
    print("Error occured while processing sample "+sample)
finally:
    result.T.to_csv(quant_support_path)