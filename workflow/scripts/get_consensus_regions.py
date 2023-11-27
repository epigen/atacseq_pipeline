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
