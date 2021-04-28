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

consensus_peaks = bedtools.BedTool(snakemake.input[1])
consensus_peaks_df = bedtools.BedTool(snakemake.input[1]).to_dataframe().set_index('name')

results=[]

for sample in annotations.index:
    peakfile=os.path.join(results_dir,"{}".format(sample),"peaks", "{}_summits.bed".format(sample))
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
        results.append(result.T)

results = [item for item in results if item is not None]
results = pd.concat(results).T
results.to_csv(output,index_label='ID')