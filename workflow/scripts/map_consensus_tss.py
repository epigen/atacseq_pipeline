#!/bin/env python

#### libraries
import pandas as pd

# map region to gene and classify if TSS
def map_region(x):
    tmp_dist = x.loc[x['homer_Distance_to_TSS'].abs().idxmin(),'homer_Distance_to_TSS']
    if (tmp_dist>=TSS_up) and (tmp_dist<=TSS_dn):
        return x.loc[x['homer_Distance_to_TSS'].abs().idxmin(),:]
    else:
        return None

#### configurations

# input
region_annotation_path = snakemake.input["region_annotation"]
consensus_counts_path = snakemake.input["consensus_counts"]

# output
tss_counts_path = snakemake.output["tss_counts"]
tss_annot_path = snakemake.output["tss_annot"]

# parameters
TSS_up = -snakemake.config["proximal_size_up"]
TSS_dn = snakemake.config["proximal_size_dn"]

# load annotations and consensus counts
annot_regions = pd.read_csv(region_annotation_path)
consensus_counts = pd.read_csv(consensus_counts_path, index_col=0)

# map regions to closest regions and subset for regions in TSS proximity
TSS_regions = annot_regions.reset_index().groupby('homer_Nearest_Ensembl').apply(map_region)
TSS_regions = TSS_regions.dropna(axis=0, how='all')

# subset the consensus counts by the successfully mapped consenesus regions, rename index to genes and save
TSS_counts = consensus_counts.loc[TSS_regions["peak_id"],:]
TSS_counts.index = TSS_regions.index
TSS_counts.to_csv(tss_counts_path)

# subset the consensus annotation by the successfully mapped consenesus regions, rename index to genes and save
annot_regions.set_index('peak_id', inplace=True)
TSS_annot = annot_regions.loc[TSS_regions["peak_id"],:]
TSS_annot.reset_index(inplace=True)
TSS_annot.index = TSS_regions.index
TSS_annot.to_csv(tss_annot_path)