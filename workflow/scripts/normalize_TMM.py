#!/bin/env python

#### normalize data with Trimmed Mean M-value Normalization (TMM) ####

#### libraries
import pybedtools
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# utility functions
def calc_norm_factors(counts_df, method):
    assert method in ['TMM', 'TMMwsp', 'RLE', 'upperquartile']
    from rpy2.robjects import numpy2ri, pandas2ri, r
    numpy2ri.activate()
    pandas2ri.activate()
    r.source(os.path.join('workflow','scripts','utils_calcNormFactors.R'))
    return r.calcNormFactors(counts_df, method=method)

def tpm(counts_matrix, transcript_len=1, norm_factors=1, log=False, pseudocount=0.5, libsize_pseudocount=1, million=1e6):
    """
    Transcript per million normalization of gene expression
    If transcript_len=1 then you get counts per million
    Setting pseudocount=0.5, libsize_pseudocount=1 ensures results equivalent to LIMMA-voom
    """
    rpk = counts_matrix.astype(float) / transcript_len
    tpm = million * (rpk + pseudocount) / (rpk.sum(axis=0) * norm_factors + libsize_pseudocount)

    return np.log2(tpm) if log else tpm

#### configurations

# input
filtered_counts=pd.read_csv(snakemake.input[0], index_col=0)

# output
output_data=snakemake.output[0]

# normalize filtered data by selected standard method
norm_method = 'TMM'

# determine norm factors
norm_factors = calc_norm_factors(filtered_counts, method=norm_method)
norm_factors = pd.Series(norm_factors, index=filtered_counts.columns, name='norm_factor')

# calculate normalized counts
normalized_lcpm = tpm(filtered_counts, norm_factors=norm_factors, log=True, pseudocount=0.5, libsize_pseudocount=1)

# save normalized counts
normalized_lcpm.to_csv(os.path.join(output_data))

# save norm factors?
# norm_factors.to_csv(os.path.join())