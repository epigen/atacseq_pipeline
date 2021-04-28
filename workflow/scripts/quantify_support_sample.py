#!/bin/env python

import os
import pandas as pd
import pybedtools as bedtools

# configuration
output = snakemake.output[0]
chrom_file = snakemake.params["chrom_file"]
sample=snakemake.wildcards["sample"]

consensus_peaks = bedtools.BedTool(snakemake.input[0])
consensus_peaks_df = bedtools.BedTool(snakemake.input[0]).to_dataframe().set_index('name')

peakfile=snakemake.input[1]
result = pd.DataFrame(0,index=consensus_peaks_df.index,columns=[sample])
    
try:
    if (peakfile is not None):
        sample_peaks = bedtools.BedTool(peakfile)
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
    result.T.to_csv(output)