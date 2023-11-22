#!/bin/env python

#### libraries
import os
import pandas as pd
import pybedtools as bedtools

#### configurations

# input
peakfiles = snakemake.input["summits_bed"]
blacklist_file = snakemake.config["blacklisted_regions"]
chrom_file = snakemake.config["chromosome_sizes"]

# output
consensus_regions_path = snakemake.output["consensus_regions"]
# results_dir = snakemake.params["results_dir"]

# parameters
slop_extension=snakemake.config["slop_extension"]


# load summits and generate consensus regions using (py)bedtools
output_bed = None

for peakfile in peakfiles:
    peak_bed = bedtools.BedTool(peakfile)
    if (blacklist_file is not None):
        peak_bed=peak_bed.intersect(blacklist_file,v=True, wa=True)

    peak_bed = peak_bed.slop(g=chrom_file, b=slop_extension)

    if (output_bed is None):
        output_bed = peak_bed
    else:
        output_bed = output_bed.cat(peak_bed,force_truncate=True)

output_bed.saveas(consensus_regions_path)
peaks = bedtools.BedTool(consensus_regions_path).sort(faidx=chrom_file).to_dataframe(names=['CHR','START','END'],dtype={'START':int,'END':int})
peaks['ID'] = peaks.index.format(formatter=(lambda x: "CONS{:011d}".format(x)))

# save results
bedtools.BedTool().from_dataframe(peaks).saveas(consensus_regions_path)


# load annotation file and get sample names
# annotations = pd.read_csv(snakemake.input[0], index_col=0)
# annotations=annotations[(annotations['pass_qc']>0)]

# # save filtered annotation file
# annotations.to_csv(snakemake.output[1])

# peakfiles = [os.path.join(results_dir,"{}".format(sample),"peaks", "{}_summits.bed".format(sample)) for sample in annotations.index]



# if not os.path.exists(os.path.split(output)[0]):
#     os.mkdir(os.path.split(output)[0])