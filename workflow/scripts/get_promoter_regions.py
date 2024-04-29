#!/bin/env python

#### libraries
import pybedtools as bedtools
import pandas as pd

# extract promoter regions
def get_promoter(feature, upstream, downstream, chrom_sizes):
    if feature.strand == '+':
        start = feature.start - upstream
        end = feature.start + downstream
    else:
        start = feature.end - downstream
        end = feature.end + upstream
    
    # Ensure start is not negative
    start = max(start, 0)
    # Ensure end is not longer than chromosome
    end = min(end, chrom_sizes.get(feature.chrom, end))
    
    # Extract gene_id and remove version number if present
    gene_id = feature.attrs['gene_id'].split('.')[0]
    
    # Create a new feature with the promoter region coordinates
    promoter = bedtools.create_interval_from_list([
        feature.chrom,
        start,
        end,
        gene_id,
        feature.attrs['gene_name'] if 'gene_name' in feature.attrs else feature.attrs['gene_id'],
#         '.',
#         feature.strand
    ])
    
    return promoter

#### configurations

# input
gtf_file = snakemake.config["gencode_gtf"]
chrom_file = snakemake.config["chromosome_sizes"]

# output
promoter_regions_path = snakemake.output["promoter_regions"]
promoter_annot_path = snakemake.output["promoter_annot"]

# parameters
TSS_up = snakemake.config["proximal_size_up"]
TSS_dn = snakemake.config["proximal_size_dn"]
genome_fasta_path = snakemake.config["genome_fasta"]

# load the genome annotation file using pybedtools
gtf = bedtools.BedTool(gtf_file)

# load and get list of valid chromosomes and sizes
chrom_sizes = {}
with open(chrom_file, 'r') as f:
    for line in f:
        chrom, size = line.strip().split('\t')
        chrom_sizes[chrom] = int(size)

# filter for features that are genes AND not pseudoautosomal regions denoted by "PAR" and create promoters
# https://www.ensembl.org/info/genome/genebuild/human_PARS.html
promoters = gtf.filter(lambda x: (x[2] == 'gene') & ("PAR" not in x["gene_id"])).each(get_promoter, TSS_up, TSS_dn, chrom_sizes)

# filter for valid chromosomes
promoters = promoters.filter(lambda x: x.chrom in chrom_sizes)

# sort promoter regions
promoters = promoters.sort(faidx=chrom_file)

# save the promoter regions to a BED file
promoters.saveas(promoter_regions_path)

# calculate GC content and length for each region and save as annotation
gc_content_length = promoters.nucleotide_content(fi=genome_fasta_path).to_dataframe()
gc_content_length.columns = [col.split('_', 1)[-1].replace('at', 'AT').replace('gc', 'GC').replace('oth', 'otherBases') for col in gc_content_length.columns]
gc_content_length = gc_content_length.add_prefix('bedtools_')
gc_content_length.columns = ["chr", "start", "end", "gene", "gene_name"] + gc_content_length.columns[5:].tolist()
gc_content_length.set_index("gene", inplace=True)
gc_content_length.to_csv(promoter_annot_path)

# load, remove last column (gene name) and save again as final promoter BED file for quantification
bedtools.BedTool(promoter_regions_path).cut(range(0, 4)).saveas(promoter_regions_path)
