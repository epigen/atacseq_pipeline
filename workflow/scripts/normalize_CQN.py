#!/bin/env python

#### normalize data with Conditional quantile normalization (CQN) ####

#### libraries
import pybedtools
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# utility functions
def bed_to_index(df):
    import pybedtools

    if isinstance(df, pybedtools.BedTool):
        df = df.to_dataframe()
    elif isinstance(df, str):
        df = pybedtools.BedTool(df).to_dataframe()
    cols = ["chrom", "start", "end"]
    if not all([x in df.columns for x in cols]):
        raise AttributeError(
            "DataFrame does not have '{}' columns.".format("', '".join(cols))
        )
    index = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(int).astype(str)
        + "-"
        + df["end"].astype(int).astype(str)
    )
    return pd.Index(index, name="region")


def get_peak_gccontent_length(bed_file, fasta_file) -> pd.DataFrame:
    import pybedtools

#     sites = pybedtools.BedTool(bed_file) if bed_file is path
    sites = pybedtools.BedTool(bed_file, from_string=True)
    nuc = sites.nucleotide_content(fi=fasta_file).to_dataframe(comment="#")[
        ["score", "blockStarts"]
    ]
    nuc.columns = ["gc_content", "length"]
    nuc.index = bed_to_index(sites)

    return nuc


def cqn(matrix, gc_content, lengths):
    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    importr("cqn")

    cqn_out = r.cqn(matrix, x=gc_content, lengths=lengths)

    y_r = cqn_out[list(cqn_out.names).index("y")]
    y = pd.DataFrame(np.array(y_r), index=matrix.index, columns=matrix.columns)
    offset_r = cqn_out[list(cqn_out.names).index("offset")]
    offset = pd.DataFrame(
        np.array(offset_r), index=matrix.index, columns=matrix.columns
    )

    return y + offset

#### configurations

# input
filtered_counts=pd.read_csv(snakemake.input[0], index_col=0)
consensus_regions=pd.read_csv(snakemake.input[1], sep='\t', index_col=3, header=None)
filtered_regions=consensus_regions.loc[filtered_counts.index,:]

# parameters
genome_fasta = snakemake.params["genome_fasta"]

# output
output_data=snakemake.output[0]

# normalize the filtered count matrices by cqn method with gc-content as covariate
# get regions
bed_regions=filtered_regions.to_string(header=False, index=False)

# get nuc info
nuc = get_peak_gccontent_length(bed_regions, fasta_file=genome_fasta)

# normalize via CQN method
normalized_counts = cqn(filtered_counts, gc_content=nuc["gc_content"], lengths=nuc["length"])

# save normalized matrix
normalized_counts.to_csv(output_data)