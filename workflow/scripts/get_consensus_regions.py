#!/bin/env python

import os
import pandas as pd
import pybedtools as bedtools

# configuration
sloop_extension=snakemake.params["sloop_extension"]
output = snakemake.output[0]
blacklist_file= snakemake.params["blacklist_file"]
chrom_file = snakemake.params["chrom_file"]
results_dir = snakemake.params["results_dir"]
# load annotation file and get sample names
annotations = pd.read_csv(snakemake.input[0], index_col=0)
annotations=annotations[(annotations['pass_qc']>0)]
peakfiles = [os.path.join(results_dir,"{}".format(sample),"peaks", "{}_summits.bed".format(sample)) for sample in annotations.index]

output_bed = None

if not os.path.exists(os.path.split(output)[0]):
    os.mkdir(os.path.split(output)[0])

for peakfile in peakfiles:
    peak_bed = bedtools.BedTool(peakfile)
    if (blacklist_file is not None):
        peak_bed=peak_bed.intersect(blacklist_file,v=True, wa=True)

    peak_bed = peak_bed.slop(g=chrom_file, b=sloop_extension)

    if (output_bed is None):
        output_bed = peak_bed
    else:
        output_bed = output_bed.cat(peak_bed,force_truncate=True)

output_bed.saveas(output)

peaks = bedtools.BedTool(output).sort(faidx=chrom_file).to_dataframe(names=['CHR','START','END'],dtype={'START':int,'END':int})
peaks['ID'] = peaks.index.format(formatter=(lambda x: "CONS{:011d}".format(x)))
bedtools.BedTool().from_dataframe(peaks).saveas(output)